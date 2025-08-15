from collections import defaultdict
import sys
import rtree
#modifed checker function to return opennets and violating nets
from checker import Inst, Net, loadAndCheck, skipCells, skipNets, markUnusedPins, adjLayer, layerColors
import LEFDEFParser
from LEFDEFParser import Rect
import numpy as np
import math
import heapq as hq
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import numpy as np
import logging

logging.basicConfig(
    filename='./src/debug_log.txt',
    level=logging.DEBUG,  # Log all messages (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(levelname)s - %(message)s'
)

layerOrient = { 'li1': 'VERTICAL', 'met1': 'HORIZONTAL', 'met2': 'VERTICAL', 'met3': 'HORIZONTAL', 'met4': 'VERTICAL', 'met5': 'HORIZONTAL' }
layerSpacing = {'li1': 170, 'met1': 140, 'met2': 140, 'met3': 300, 'met4': 300, 'met5': 1600}
layerWidth = {'li1': 170, 'met1': 140, 'met2': 140, 'met3': 300, 'met4': 300, 'met5': 1600}



def plot_net_all_layers(net_name, nets, guide_data, track_data=None, via_points_data=None, vertices=None, blocked_vertices = None, new_sources = None):
    target_net = next((net for net in nets if net._name == net_name), None)
    if not target_net:
        print(f"Net {net_name} not found.")
        return

    fig, ax = plt.subplots(figsize=(10, 10))

    # Assign distinct colors per layer
    all_layers = sorted(track_data.keys()) if track_data else []
    layer_cmap = cm.get_cmap('tab20', len(all_layers))
    layer_to_color = {layer: layer_cmap(i) for i, layer in enumerate(all_layers)}

    # ----- Draw tracks in different faint colors -----
    if track_data:
        for layer, lines in track_data.items():
            color = layer_to_color[layer]
            for line in lines:
                x_coords, y_coords = zip(*line)
                ax.plot(x_coords, y_coords, color=color, alpha=0.8, linewidth=0.7, label=f'Track: {layer}')

    # ----- Draw pin shapes -----
    for pin_name, shapes in target_net._pins.items():
        for layer, rects in shapes.items():
            for r in rects:
                x1, y1 = r.ll.x, r.ll.y
                x2, y2 = r.ur.x, r.ur.y
                width = abs(x2 - x1)
                height = abs(y2 - y1)
                ll = (min(x1, x2), min(y1, y2))
                color = layer_to_color.get(layer, 'blue')
                rect_patch = patches.Rectangle(
                    ll, width, height,
                    edgecolor=color, facecolor=color, alpha=0.3,
                    label=f'Pin: {layer}'
                )
                ax.add_patch(rect_patch)

    # ----- Draw guide rectangles -----
    if net_name in guide_data:
        for x_ll, y_ll, x_ur, y_ur, layer in guide_data[net_name]:
            width = x_ur - x_ll
            height = y_ur - y_ll
            color = layer_to_color.get(layer, 'red')
            guide_patch = patches.Rectangle(
                (x_ll, y_ll), width, height,
                edgecolor=color, facecolor=color, alpha=0.2,
                label=f'Guide: {layer}'
            )
            ax.add_patch(guide_patch)

    # ----- Plot all escape via points -----
    if via_points_data and net_name in via_points_data:
        plotted = set()
        for pin_name, layer_map in via_points_data[net_name].items():
            for src_layer, pts in layer_map.items():
                for x, y, tgt_layer in pts:
                    label = f'Via to {tgt_layer}' if tgt_layer not in plotted else ''
                    ax.plot(x, y, 'o', markersize=4,
                            markerfacecolor=layer_to_color.get(tgt_layer, 'black'),
                            markeredgecolor='k', label=label)
                    plotted.add(tgt_layer)

    # draw all vertices in the guide
    
        for x, y,layer in vertices:
                       # plot the vertices
            ax.plot(x, y, 'x', markersize=4, color='black', label='Guide Vertices')
       
    #draw blocked vertices
    if blocked_vertices:
        for x,y,layer in blocked_vertices:
            ax.plot(x,y,'+',markersize = 5,color = 'red',label='blocked_vertex')
            
    #new source:
    if new_sources:
        for x,y,layer in new_sources:
            ax.plot(x,y,'o',markersize=5,color='green',label='new_sources')
    # ----- Final Formatting -----
    ax.set_aspect('equal')
    ax.set_title(f"All Layers View for Net {net_name}")
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")

    # Unique legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.show()

def verify_net_rects_in_rtree(net, rtree_cache, layerwidth, layerspacing):
    """
    Checks if all bloated rectangles in net._sol exist in the R-tree.
    """
    missing = []
    for layer, rects in net._sol.items():
        if layer not in rtree_cache:
            print(f"[!] No R-tree found for layer {layer}")
            continue

        tree = rtree_cache[layer]
        spacing = layerspacing[layer] + (layerwidth[layer] // 2)

        for r in rects:
            bbox = (r.ll.x - spacing, r.ll.y - spacing, r.ur.x + spacing, r.ur.y + spacing)
            hits = list(tree.intersection(bbox, objects=True))

            found = any(
                obj.object.get("net") == net._name
                for obj in hits if isinstance(obj.object, dict)
            )

            if not found:
                missing.append((layer, bbox))

    if not missing:
        print(f"[✓] All rectangles of net {net._name} are present in the R-tree.")
    else:
        print(f"[❌] Missing {len(missing)} rects from R-tree for net {net._name}:")
        for layer, bbox in missing:
            print(f"  -> Layer {layer}, bbox: {bbox}")
class GUIDEParser:
    def __init__(self, guide_file):
        self.guide_file = guide_file
        self.nets = {}
        self.parse()
        

    def parse(self):
        
        with open(self.guide_file, 'r') as file:
            current_net = None
            inside_paraenthesis = False
            for lines in file:
                line = lines.strip()
                if not line or line.startswith("#"):
                    continue
                if line == "(":
                    inside_paraenthesis = True
                elif line == ")":
                    inside_paraenthesis = False
                elif inside_paraenthesis and current_net:
                    parts = line.split()
                    if len(parts) == 5:
                        x_ll, y_ll, x_ur, y_ur, layer = parts
                        x_ll, y_ll, x_ur, y_ur = map(int, (x_ll, y_ll, x_ur, y_ur))
                        segment = (x_ll, y_ll, x_ur, y_ur, layer)
                        self.nets[current_net].append(segment)
                elif not inside_paraenthesis and not line.startswith("(") and not line.startswith(")"):
                    # Line is a net name
                    current_net = line
                    self.nets[current_net] = []
        return self.nets
class Vertex:
  def __init__(self, x, y, l, cost=math.inf, parent=None, nbrs=None):
    self._xyl = (x, y,l)
    self._cost = cost
    self._parent = parent
    self._nbrs = nbrs
  def __lt__(self, r):
    return self._cost < r._cost
  def __eq__(self, r):
    return self._xyl == r._xyl
  def __repr__(self):
    return f'({self._xyl})'
class priority_queue:
  def __init__(self, vertices = []):
    self._vertices = vertices[:]
    self._q = vertices[:]
    hq.heapify(self._q)
  def push(self, v):
    hq.heappush(self._q, v)
  def pop(self):
    return(hq.heappop(self._q))
  def update(self, v, cost):
    try: i = self._q.index(v)
    except ValueError: i = None
    if i is not None:
      self._q[i]._cost = cost
      hq.heapify(self._q)
  def updateIndex(self, i, cost):
    assert i < len(self._q)
    self._vertices[i]._cost = cost
    hq.heapify(self._q)
  def empty(self):
    return len(self._q) == 0
  def __contains__(self, v):
    return v in self._q
  def __repr__(self):
    return str(self._q)
def dist(u, v):
  via_cost = 10
  dz = via_cost if u._xyl[2] != v._xyl[2] else 0
  return abs(u._xyl[0] - v._xyl[0]) + abs(u._xyl[1] - v._xyl[1]) + dz

def astar(V, s, t):
    #print("starting astar")
    for v in V:
        v._cost, v._parent = math.inf, None
    
    vertex_dict = {id(v): v for v in V}

    # Dictionary to track g_values and h_values
    g_values = {id(v): math.inf for v in V}
    h_values = {id(v): dist(v, t) for v in V}
    g_values[id(s)] = 0

    # Set f_value for start vertex (g+h)
    s._cost = g_values[id(s)] + h_values[id(s)]

    #Tie breaking: when f_values are equal, higher g_values gets priority
    #for this substract a small epsilon times g_value from the overall cost

    epsilon = 1e-10
    s._cost -= epsilon * g_values[id(s)]
    # Initialize priority queue
    Q = priority_queue([s])

    while not Q.empty():
        u = Q.pop()
        #print(f"expanding: {u}  with cost: {u._cost} ")

        # If target reached, break
        if u == t:
            #print("reached target: {u}")
            break

        # Explore neighbors
        for v in u._nbrs:            
                new_g = g_values[id(u)] + dist(u, v)
                if new_g < g_values[id(v)]:
                    v._parent = u
                    g_values[id(v)] = new_g
                    v._cost = (new_g + h_values[id(v)]) - epsilon * new_g
                    # If v is in Q, update its cost
                    if v not in Q:
                        Q.push(v)
                    else:
                        Q.update(v, v._cost)
    if not t._parent:
        #print("No valid path found")
        return None

    # Reconstruct path
    path = [t]
    while path[-1]._parent is not None:
        path.append(path[-1]._parent)
    return path

def decompose_net(net):
    """Return list of (pin1, pin2) pairs using Prim's algorithm"""
    pins = list(net._pins.items())  # [(pin_name, shapes), ...]
    
    if len(pins) < 2:
        return []

    # Build complete distance graph
    graph = defaultdict(dict)
    for i, (name1, shapes1) in enumerate(pins):
        for j, (name2, shapes2) in enumerate(pins[i+1:]):
            # Find closest points between any layers
            min_dist = float('inf')
            for layer1, rects1 in shapes1.items():
                for layer2, rects2 in shapes2.items():
                    c1 = (rects1[0].xcenter(), rects1[0].ycenter())  # Fixed here
                    c2 = (rects2[0].xcenter(), rects2[0].ycenter())  # And here
                    
                    # Add layer transition cost
                    dist = abs(c1[0]-c2[0]) + abs(c1[1]-c2[1]) 
                    if layer1 != layer2:
                        dist += layerSpacing[layer1]* 3
                        
                    min_dist = min(min_dist, dist)
                    
            graph[name1][name2] = min_dist
            graph[name2][name1] = min_dist


    # Prim's algorithm for MST
    mst = []
    visited = set([pins[0][0]])
    while len(visited) < len(pins):
        min_edge = (None, None, float('inf'))
        for u in visited:
            for v, weight in graph[u].items():
                if v not in visited and weight < min_edge[2]:
                    min_edge = (u, v, weight)
        if min_edge[0] is None:
            break
        mst.append((min_edge[0], min_edge[1]))
        visited.add(min_edge[1])
        
    return mst

def calculate_hpwl(net,subnets):
    # all_points = []
    # for pin in net._pins.values():
    #     for layer, rects in pin.items():
    #         for rect in rects:
    #             # Consider layer changes in HPWL
    #             all_points.append((
    #                 rect.xcenter(), 
    #                 rect.ycenter(),
    #                 layer  # Include layer in HPWL calculation
    #             ))
    
    # if not all_points:
    #     return 0
    
    # # Find bounding box across all layers
    # min_x = min(p[0] for p in all_points)
    # max_x = max(p[0] for p in all_points)
    # min_y = min(p[1] for p in all_points)
    # max_y = max(p[1] for p in all_points)
    
    # return (max_x - min_x) + (max_y - min_y)
    
    total_hpwl = 0
    for pin1, pin2 in subnets[net._name]:
        pin1_shapes = net._pins.get(pin1, {})
        pin2_shapes = net._pins.get(pin2, {})

        pin1_coords = []
        pin2_coords = []

        for rects in pin1_shapes.values():
            for r in rects:
                pin1_coords.append((r.xcenter(), r.ycenter()))

        for rects in pin2_shapes.values():
            for r in rects:
                pin2_coords.append((r.xcenter(), r.ycenter()))

        if not pin1_coords or not pin2_coords:
            continue

        # Use first center point from each pin (or compute average)
        x1, y1 = pin1_coords[0]
        x2, y2 = pin2_coords[0]

        hpwl = abs(x1 - x2) + abs(y1 - y2)
        total_hpwl += hpwl

    return total_hpwl 

def check_obstacles_in_guide(net_name, guide_data, rtree_cache, ignore_self=True):
    print(f"\n[Guide Check] Net: {net_name}")
    
    guide_rects_by_layer = {}
    for (x1, y1, x2, y2, layer) in guide_data.get(net_name, []):
        r = Rect(x1, y1, x2, y2)
        guide_rects_by_layer.setdefault(layer, []).append(r)

    total_hits = 0
    for layer, rects in guide_rects_by_layer.items():
        if layer not in rtree_cache:
            print(f"  [!] No R-tree for layer {layer}")
            continue

        for r in rects:
            bbox = (r.ll.x, r.ll.y, r.ur.x, r.ur.y)
            hits = list(rtree_cache[layer].intersection(bbox, objects=True))
            true_hits = []

            for item in hits:
                obj = item.object
                if not isinstance(obj, dict):
                    continue
                if ignore_self and obj.get("net") == net_name:
                    continue
                true_hits.append(obj)

            if true_hits:
                print(f"  [X] Obstacles in {bbox} on {layer}:")
                for obj in true_hits:
                    print(f"     - {obj}")
                total_hits += len(true_hits)
            else:
                print(f"  [✓] No obstacles in {bbox} on {layer}")

    print(f"  --> Total obstacle hits: {total_hits}")
 
def get_escape_via_points(rect, layer, track_data, layer_orient, adj_layer_map):
    """
    Return list of escape via points for a given shape on a layer.
    For li1: elevate to met1.
    For met1 to met5: use intersection points with adjacent layer tracks inside the shape bounds.
    """
    via_points = []

    x1, y1 = rect.ll.x, rect.ll.y
    x2, y2 = rect.ur.x, rect.ur.y
    x_min, x_max = sorted([x1, x2])
    y_min, y_max = sorted([y1, y2])

    if layer == "li1":
        # Elevate to met1 by overlapping li1 (VERTICAL) and met1 (HORIZONTAL) tracks
        li1_tracks = track_data.get("li1", [])
        met1_tracks = track_data.get("met1", [])

        if not li1_tracks or not met1_tracks:
            return []

        li1_xs = [line[0][0] for line in li1_tracks if x_min <= line[0][0] <= x_max]
        met1_ys = [line[0][1] for line in met1_tracks if y_min <= line[0][1] <= y_max]

        for x in li1_xs:
            for y in met1_ys:
                if x_min <= x <= x_max and y_min <= y <= y_max:
                    via_points.append((x, y, "met1"))
        return via_points

    # For routing layers: check intersections with adjacent layers
    if layer not in track_data or layer not in layer_orient:
        return []

    base_tracks = track_data[layer]
    base_orient = layer_orient[layer]

    for adj in adj_layer_map.get(layer, []):
        if adj not in track_data or adj not in layer_orient:
            continue
        adj_tracks = track_data[adj]
        adj_orient = layer_orient[adj]

        if base_orient == 'HORIZONTAL' and adj_orient == 'VERTICAL':
            base_ys = [line[0][1] for line in base_tracks if y_min <= line[0][1] <= y_max]
            adj_xs = [line[0][0] for line in adj_tracks if x_min <= line[0][0] <= x_max]

            for y in base_ys:
                for x in adj_xs:
                    if x_min <= x <= x_max and y_min <= y <= y_max:
                        via_points.append((x, y, adj))

        elif base_orient == 'VERTICAL' and adj_orient == 'HORIZONTAL':
            base_xs = [line[0][0] for line in base_tracks if x_min <= line[0][0] <= x_max]
            adj_ys = [line[0][1] for line in adj_tracks if y_min <= line[0][1] <= y_max]

            for x in base_xs:
                for y in adj_ys:
                    if x_min <= x <= x_max and y_min <= y <= y_max:
                        via_points.append((x, y, adj))

    return via_points

def extract_vertices_from_guide(net_name, guide_data, track_data, layer_orient, expand = False,iter=1):
    #vertices = {}
    vertices= set()
    adjLayer = {'met1': ['met2'], 'met2': ['met3'], 'met3': ['met4'], 'met4': ['met5']}
    guide_rects = guide_data[net_name]
    
    for x1, y1, x2, y2, layer in guide_rects:
        if expand:
            inflate = iter*1000 #4000 microns
            x1,y1,x2,y2 = x1-inflate, y1-inflate, x2 +inflate, y2 +inflate
        
        if layer == 'li1':
            continue

        for adj in adjLayer.get(layer, []):
            
            # has_guide_on_adj = any(g_layer == adj for *_, g_layer in guide_data[net_name])
            # if not has_guide_on_adj:
            #         continue
            
            # if adj == 'li1':
            #     continue

            if layer not in track_data or adj not in track_data:
                continue

            # Determine which is H and which is V
            if layer_orient[layer] == 'HORIZONTAL':
                h_layer = layer
                v_layer = adj
            else:
                h_layer = adj
                v_layer = layer

            # Get horizontal and vertical track coords
            h_tracks = [line[0][1] for line in track_data[h_layer] if y1 <= line[0][1] <= y2]
            v_tracks = [line[0][0] for line in track_data[v_layer] if x1 <= line[0][0] <= x2]

            for x in v_tracks:
                for y in h_tracks:
                    
                    vertices.add((x, y, layer))
                    vertices.add((x, y, adj))
    

    return vertices

def filter_vertices_by_rtree(vertices, rtree_cache, current_net_name):
    
    layerwidth =  {'li1': 170, 'met1': 140, 'met2': 140, 'met3': 300, 'met4': 300, 'met5': 1600}
    layerspacing =  {'li1': 170, 'met1': 140, 'met2': 140, 'met3': 300, 'met4': 300, 'met5': 1600}
    blocked = set()
    for x, y, layer in vertices:
        if layer not in rtree_cache:
            continue
        
        
        epsilon = layerspacing[layer] +layerwidth[layer]//2
        #epsilon = 1
        box = (x-epsilon, y-epsilon, x+epsilon, y+epsilon)
        for item in rtree_cache[layer].intersection(box, objects=True):
            obj = item.object
            if isinstance(obj, dict):
                if obj.get('type') == 'net' and obj.get('net') == current_net_name:
                    continue  # Ignore this net's own shapes
                #print(f"blocked at ({x},{y},{layer}) by {obj}")
            blocked.add((x, y, layer))
            
                
            break  # no need to check further
            
    return blocked

def build_3d_grid_graph(vertices, layer_orient, escape_points=None):
    import networkx as nx
    from collections import defaultdict

    G = nx.Graph()
    vertex_map = defaultdict(list)  # key: layer -> list of (x,y) in that layer

    # Organize vertices by layer for easier access
    for x, y, layer in vertices:
        vertex_map[layer].append((x, y))
        G.add_node((x, y, layer))

    # Intra-layer edges
    for layer, points in vertex_map.items():
        orient = layer_orient[layer]
        if orient == 'HORIZONTAL':
            # Same Y, sorted X
            line_groups = defaultdict(list)
            for x, y in points:
                line_groups[y].append(x)
            for y, x_vals in line_groups.items():
                x_vals.sort()
                for i in range(len(x_vals) - 1):
                    u = (x_vals[i], y, layer)
                    v = (x_vals[i + 1], y, layer)
                    G.add_edge(u, v)
        elif orient == 'VERTICAL':
            # Same X, sorted Y
            line_groups = defaultdict(list)
            for x, y in points:
                line_groups[x].append(y)
            for x, y_vals in line_groups.items():
                y_vals.sort()
                for i in range(len(y_vals) - 1):
                    u = (x, y_vals[i], layer)
                    v = (x, y_vals[i + 1], layer)
                    G.add_edge(u, v)

    # Inter-layer (via) edges
    for x, y, layer in vertices:
        u = (x, y, layer)
        for adj in adjLayer.get(layer, []):
            v = (x, y, adj)
            if v in G:
                G.add_edge(u, v)
                G.add_edge(v, u)  #bidirectional edg

    # Connect escape points (if provided)
    if escape_points:
        for net_pin, point_dict in escape_points.items():
            for layer, points in point_dict.items():
                for x, y, _ in points:
                    if (x, y, layer) not in G:
                        G.add_node((x, y, layer))
                    for nbr in vertex_map[layer]:
                        if nbr == (x, y):
                            G.add_edge((x, y, layer), (x, y, layer))  # trivial connection or replace with routing rule

    return G

def expand_sources_with_path(path,vertex_set):
   if not path:
       return
   #print("expanding source list")
   #stride = max(1,len(path)//max_count)
   #for i in range(0,len(path),stride):
   for i in range(0,len(path)):
       vertex_set.add(path[i]._xyl)

def buildTree(nets, insts, obsts):
  import rtree
  lT = {layer: rtree.index.Index() for layer in layerColors}
  obstid = len(nets)

  count = 0
  for inst in insts.values():
    for layer, rects in inst._obsts.items():
      for r in rects:
        lT[layer].insert(count, (r.ll.x, r.ll.y, r.ur.x, r.ur.y), obj=obstid)
        count += 1

  for layer, rects in obsts.items():
    for r in rects:
      lT[layer].insert(count, (r.ll.x, r.ll.y, r.ur.x, r.ur.y), obj=obstid)
      count += 1

  for net in nets:
    for layer, rects in net._sol.items():
      for r in rects:
        lT[layer].insert(count, (r.ll.x, r.ll.y, r.ur.x, r.ur.y), obj=net._id)
        count += 1

    for p, lr in net._pins.items():
      for layer, rects in lr.items():
        for r in rects:
          lT[layer].insert(count, (r.ll.x, r.ll.y, r.ur.x, r.ur.y), obj=net._id)
          count += 1

  return lT

def convert_path_to_rects(path, layer_width):
    from LEFDEFParser import Rect
    segments_by_layer = {}

    i = 0
    while i < len(path) - 1:
        
        x1, y1, l1 = path[i]._xyl
        x2, y2, l2 = path[i + 1]._xyl

        # If it's a via (layer change at the same x,y), insert via box only at then end and begiing of path
        if l1 != l2 and x1 == x2 and y1 == y2:
            if(i == 0 or i == len(path)-2):
                for l in [l1, l2]:
                    w = layer_width[l] // 2 
                    #s = layer_spacing[l]
                    r = Rect(x1 - w, y1 - w, x1 + w, y1 + w)
                    if l not in segments_by_layer:
                        segments_by_layer[l] = []
                    segments_by_layer[l].append(r)
            i += 1
            continue
          
        # Collapse collinear segments on same layer
        j = i + 1
        while j < len(path) - 1:
            xn, yn, ln = path[j + 1]._xyl
            if ln != l1:
                break
            if not ((x2 == xn and y2 != yn) or (y2 == yn and x2 != xn)):
                break
            x2, y2 = xn, yn
            j += 1

        # Bloat by half-width
        w = layer_width[l1] // 2
        x_min = min(x1, x2) - w
        x_max = max(x1, x2) + w
        y_min = min(y1, y2) - w
        y_max = max(y1, y2) + w
        r = Rect(x_min, y_min, x_max, y_max)
        if l1 not in segments_by_layer:
            segments_by_layer[l1] = []
        segments_by_layer[l1].append(r)
        i = j

    return segments_by_layer

def add_rects_to_rtree(tree, net_id, net_sol, layerwidth, layerspacing, start_id):
    
    count = start_id
    #bloated_rects = defaultdict(list)
    
    for layer, rects in net_sol.items():
        s = (layerwidth[layer]//2) + layerspacing[layer]
        
        for r in rects:
            tree[layer].insert(count,(r.ll.x-s, r.ll.y-s, r.ur.x + s, r.ur.y + s),obj=net_id)
            count += 1
            
            
    # for path in path_list:
    #     segments = convert_path_to_rects(path[1:-1], layerwidth)
    #     for layer, rects in segments.items():
    #         bloated_rects[layer].extend(rects)

    # for layer, rects in bloated_rects.items():
    #     if layer not in tree:
    #         tree[layer] = rtree.index.Index()
    #     spacing = layerspacing[layer] + (layerwidth[layer]//2)
    #     for r in rects:
    #         #print(layer)
    #         tree[layer].insert(count, (r.ll.x - spacing, r.ll.y - spacing, r.ur.x + spacing, r.ur.y + spacing),  obj={'type': 'net', 'net': net_id})
    #         count += 1

    return count

def write_paths_to_def(d, output_path, all_net_paths, layer_width):
    
    #print(f"inside")
    #print(len(all_net_paths))
    for net in d.nets():
        if net.name() not in all_net_paths:
            #print("nonet")
            continue
        for path in all_net_paths[net.name()]:
            #print(f"path to be converterd: {path}")
            segments = convert_path_to_rects(path[1:-1], layer_width)  # remove dummy source/target
            #print(f"segments: {segments}")
            for layer, rects in segments.items():
                for r in rects:
                    net.addRect(layer, r.ll.x, r.ll.y, r.ur.x, r.ur.y)

    d.writeDEF(output_path)

def elevate_and_snap_point(x, y, current_layer, track_data, layer_orient, adjLayer, blocked):
    """
    Attempt to elevate a point to a higher layer and snap to nearby track.
    Return the elevated (x, y, new_layer) if successful, else None.
    """
    visited = set()
    layers_to_try = adjLayer.get(current_layer, [])

    for l in layers_to_try:
        if l in visited or l not in track_data:
            continue
        visited.add(l)
        orient = layer_orient[l]

        if orient == 'HORIZONTAL':
            ys = [line[0][1] for line in track_data[l]]
            closest_y = min(ys, key=lambda y1: abs(y1 - y))
            candidate = (x, closest_y, l)
        else:
            xs = [line[0][0] for line in track_data[l]]
            closest_x = min(xs, key=lambda x1: abs(x1 - x))
            candidate = (closest_x, y, l)

        if candidate not in blocked:
            return candidate

    return None

def force_route(net_name, vertices_set, blocked, sources, targets, track, layerOrient, adjLayer, guide_data, via_points, attempt_expand=False):
    #print(f"Starting force-routing for net {net_name}...")

    def get_unblocked(vertices, allow_extra=[]):
        return [v for v in vertices if v not in blocked or v in allow_extra]

    # Case 1: Only sources are blocked
    blocked_sources = [s for s in sources if s in blocked]
    blocked_targets = [t for t in targets if t in blocked]

    unblocked_vertices = []

    if blocked_sources and not blocked_targets:
        #print(" → Unblocking sources only")
        unblocked_vertices = get_unblocked(vertices_set, allow_extra=sources)
    elif blocked_targets and not blocked_sources:
        #print(" → Unblocking targets only")
        unblocked_vertices = get_unblocked(vertices_set, allow_extra=targets)
    elif blocked_sources or blocked_targets:
        #print(" → Unblocking both sources and targets")
        unblocked_vertices = get_unblocked(vertices_set, allow_extra=sources + targets)
    else:
        # All escape points are already unblocked, something else failed
        #print(" → All sources/targets are unblocked, trying again...")
        unblocked_vertices = get_unblocked(vertices_set)

    # Build graph and try routing
    G = build_3d_grid_graph(unblocked_vertices, layerOrient, via_points[net_name])
    path = try_astar_with_virtuals(G, sources, targets)

    if path:
        #print("Force routing succeeded (with unblocking)")
        return path

    #Try elevation
    #print(" → Unblocking failed. Trying to expand grid  ...")
    if not attempt_expand:
        #print(" → Expanding guide and retrying...")
        vertices_set_exp = extract_vertices_from_guide(net_name, guide_data, track, layerOrient, expand=True, iter=3)
        path = force_route(net_name, vertices_set_exp, blocked, sources, targets, track, layerOrient, adjLayer, guide_data, via_points, attempt_expand=True)
        if path:
            #print(" orce routing succeeded (with elevation)")
            return path
    #print(" → Unblocking failed. Trying to elevate points ...")
    def elevate(points):
        elevated = []
        for x, y, layer in points:
            elev = elevate_and_snap_point(x, y, layer, track, layerOrient, adjLayer, blocked)
            if elev:
                elevated.append(elev)
        return elevated

    elevated_sources = elevate(sources)
    elevated_targets = elevate(targets)
    if not elevated_sources and not elevated_targets:
        print(" → Elevation failed. No higher track access.")
        #return []
    else:
        combined = list(vertices_set) + elevated_sources + elevated_targets
        G = build_3d_grid_graph(combined, layerOrient, via_points[net_name])
        path = try_astar_with_virtuals(G, elevated_sources or sources, elevated_targets or targets)
        if path:
            #print("Force routing succeeded (with elevation)")
            return path

    
    

    #print(" All force-routing attempts failed.")
    
    #sys.exit(1)
    return None

def try_astar_with_virtuals(G, sources, targets):
    G.add_node("V_SOURCE")
    G.add_node("V_TARGET")
    for s in sources:
        if s in G:
            G.add_edge("V_SOURCE", s, weight=0)
    for t in targets:
        if t in G:
            G.add_edge("V_TARGET", t, weight=0)

    vertex_map = {}
    for node in G.nodes():
        if node in ("V_SOURCE", "V_TARGET"):
            continue
        x, y, l = node
        v = Vertex(x, y, l)
        v._nbrs = []
        vertex_map[node] = v
    for u, v, _ in G.edges(data=True):
        if u in vertex_map and v in vertex_map:
            vertex_map[u]._nbrs.append(vertex_map[v])
            vertex_map[v]._nbrs.append(vertex_map[u])
    v_source = Vertex(-1, -1, -1, cost=0)
    v_target = Vertex(-2, -2, -2)
    v_source._nbrs = [vertex_map[n] for n in G.neighbors("V_SOURCE") if n in vertex_map]
    v_target._nbrs = [vertex_map[n] for n in G.neighbors("V_TARGET") if n in vertex_map]
    for n in G.neighbors("V_TARGET"):
        if n in vertex_map:
            vertex_map[n]._nbrs.append(v_target)
    return astar(list(vertex_map.values()) + [v_source, v_target], v_source, v_target)

def precompute_net_vertices(net_name, guide_data, track, layer_orient, rtree_cache):
    """Precompute vertices and filtered sets for a net"""
    vertices_set = extract_vertices_from_guide(net_name, guide_data, track, layer_orient)
    
    return vertices_set

def spacing_net(spacingnets):
    # Frequency and violation count per net
   
    violations = defaultdict(int)

    for a, b in spacingnets:
        if a != 'obst':
            
            violations[a] += 1
        if b != 'obst':
           
            violations[b] += 1
    

    # Nets to rip
    nets_to_rip = set()

    for a, b in spacingnets:
        if a == 'obst':
            nets_to_rip.add(b)
        elif b == 'obst':
            nets_to_rip.add(a)
        else:
            # Step 1: Violation count priority
            if violations[a] > violations[b]:
                nets_to_rip.add(a)            
            else:                
                nets_to_rip.add(b)
    #print(f"violations: {violations}")
    #print(f"Final nets selected for rip-up: {nets_to_rip}")
    return nets_to_rip

def route_net(iter,net, rtree_cache,net_subnets, guide_data, layerOrient, via_points):
    """iteration,
                net,
                rtree_cache,
                net_subnets,                
                guide_data,                
                layerOrient,
                via_points
    Attempt to reroute a net. Return success status and updated obstacle_id.
    """
    all_paths = [] #to store the path of net, list of list in case of subnets
    net_name = net._name
    subnets = net_subnets[net_name]
    if not subnets:
        return False,[]

    # Expand guide area slightly for reroute recovery    
    vertices_set = extract_vertices_from_guide(net_name, guide_data, track, layerOrient, expand=True,iter=iter+1)
    #print(f"vertices: {len(vertices_set)}")
    blocked = filter_vertices_by_rtree(vertices_set, rtree_cache,net_name)
    #print(f"blocked:{len(blocked)}")
    routed_vertices_so_far = set()

    for pin1, pin2 in subnets:
        #print(f"{pin1}->{pin2}")
        pin1_vias = via_points[net_name].get(pin1, {})
        pin2_vias = via_points[net_name].get(pin2, {})
        sources = [(x, y, _) for l, pts in pin1_vias.items() for (x, y, _) in pts]
        targets = [(x, y, _) for l, pts in pin2_vias.items() for (x, y, _) in pts]
        sources += list(routed_vertices_so_far)
        #print(f"sources:{len(sources)},targets:{len(targets)}")

        if not sources or not targets:
            print(f"Missing escape points")
            return False,[]

        unblocked_vertices = [v for v in vertices_set if v not in blocked]
        G = build_3d_grid_graph(unblocked_vertices, layerOrient, via_points[net_name])   
        path = try_astar_with_virtuals(G, sources, targets)        
        if not path:
            print(f"path not found.") 
            return False,[]                     
        else:
            print(f"path found for {pin1}->{pin2}")
            all_paths.append(path[::-1])
            print(f"path len: {len(path)}")
            expand_sources_with_path(path[1:-1], routed_vertices_so_far)

    print(f"net in reroute net")          
    
                
    

    return True,all_paths

def rrr(iter,open_nets, violating_nets,netDict,net_subnets,all_net_paths,insts, obsts, nets, guide_data, layerOrient, via_points, layerwidth, output_def, original_def, lef):
    import time
    # from utilits import convert_path_to_rects, add_rects_to_rtree
    # from checker import check

    max_iter = iter    
    max_time = 3600  # 60 minutes
    start_time = time.time()    
    #print(f"violating_nets: {violating_nets}")
    print(f"open_nets: {open_nets}")
    #flat_violating_nets = list({v[0] for v in violating_nets})  # set to remove duplicates
    print(f"violating_nets: {violating_nets}")
    #reroute_nets = set(open_nets) | set(violating_nets)
    
    for iteration in range(1, max_iter + 1):
        print(f"\n=== RRR Iteration {iteration} ===")
        start_time=time.time()
        logging.info(f"Iteration:{iteration},start_time:{start_time}")
        if time.time() - start_time > max_time:
            print("Timeout: Routing exceeded 30 minutes")
            logging.warning(f"Timeout")
            break
        
        reroute_nets = set(open_nets) | set(violating_nets)
        
        if not set(open_nets) or not set(violating_nets):
            print("No nets to reroute. Stopping early.")
            logging.info(f"Stopped early")
            break

        # Clear old routes for reroute nets
        for net in reroute_nets:           
            
                print(f"clearing initial_sol: {len(netDict[net]._sol)}")
                
                netDict[net]._sol.clear()
                print(f"clearing initial_sol: {len(netDict[net]._sol)}")

        # Rebuild R-tree with only valid nets        
        rtree_cache = buildTree(nets, insts, obsts)        
        for net in nets: #net is an object
            if net._name not in reroute_nets:
              continue  
            net_name = net._name
            #routing each net
            
            success,all_paths= route_net(
                iteration,
                net,
                rtree_cache,
                net_subnets,                
                guide_data,                
                layerOrient,
                via_points
                
            )
            
            if success:
                print("pass, path converted to rect and added to net._sol")
                all_net_paths[net_name] = all_paths  
                for path in all_paths:
                    rects = convert_path_to_rects(path[1:-1],layerWidth)
                    for layer, rs in rects.items():
                        if layer not in net._sol:
                            net._sol[layer] = []
                        net._sol[layer].extend(rs)                            
                
            else:
                print(f"failed")
        
        
        ideff = LEFDEFParser.DEFReader()
        ideff.readDEF(original_def)
        # Write updated paths to DEF
        write_paths_to_def(
            d=ideff,
            output_path=output_def,
            all_net_paths=all_net_paths,
            layer_width=layerwidth
        )

        # Rerun checker
        logging.info(f"End of iteration{iteration} time:{time.time()-start_time}")
        print(f"\nRunning DRC after RRR pass{iteration}...")
        logging.info(f"Running drc checker:after{iteration}")
        netDict,insts,pins,obsts,nets,spacingnets,openset = loadAndCheck(output_def, original_def, lef, plot=False)
        violating_nets = spacing_net(spacingnets)
        if not violating_nets and not openset:
            print("✅ All nets routed and DRC clean")
            logging.info(f"All nets routed and drc clean")
            break
        else:
            print(f"⚠️ Open nets: {len(openset)} | DRC-violating nets: {len(violating_nets)}")
            open_nets = openset
            violating_nets = violating_nets
            logging.info(f"open sets: {len(open_nets)}")

def load_design_files(lef, idef, guide):
    from LEFDEFParser import LEFReader, DEFReader
    from GUIDEParser import GUIDEParser
    from collections import defaultdict
    
    leff = LEFReader()
    leff.readLEF(lef)

    ideff = DEFReader()
    ideff.readDEF(idef)

    guidefile = GUIDEParser(guide)
    guide_data = guidefile.nets

    lefDict = {m.name(): m for m in leff.macros()}
    insts = {inst.name(): Inst(inst, lefDict[inst.macro()]) for inst in ideff.components() if inst.macro() not in skipCells}

    pins = defaultdict(lambda: defaultdict(list))
    for p in ideff.pins():
        for port in p.ports():
            for layer, rects in port.items():
                pins[p.name()][layer].extend(Rect(r.ll.x, r.ll.y, r.ur.x, r.ur.y) for r in rects)

    nets = []
    idx = 0
    for net in ideff.nets():
        if net.name() not in skipNets:
            nets.append(Net(net, insts, pins, idx))
            idx += 1

    return leff, ideff, guide_data, insts, pins, nets

def preprocess_routing_data(ideff, nets, insts, pins, guide_data):
    from collections import defaultdict

    bbox = ideff.bbox()
    track = {}
    for l, orient in layerOrient.items():
        ltracks = ideff.tracks().get(l, [])
        matched = next((t for t in ltracks if (orient == 'VERTICAL' and t.orient == 'X') or (orient == 'HORIZONTAL' and t.orient == 'Y')), None)
        if matched:
            if orient == 'VERTICAL':
                track[l] = [[(matched.x + i * matched.step, bbox.ll.y), (matched.x + i * matched.step, bbox.ur.y)] for i in range(matched.num)]
            else:
                track[l] = [[(bbox.ll.x, matched.x + i * matched.step), (bbox.ur.x, matched.x + i * matched.step)] for i in range(matched.num)]

    net_subnets = {net._name: decompose_net(net) for net in nets}
    nets.sort(key=lambda n: calculate_hpwl(n, net_subnets), reverse = False)
    netDict = {net._name: net for net in nets}

    via_points = {}
    for net in nets:
        net_name = net._name
        via_points[net_name] = {}
        for pin_name, shapes in net._pins.items():
            via_points[net_name][pin_name] = {}
            for layer, rects in shapes.items():
                pts = [pt for rect in rects for pt in get_escape_via_points(rect, layer, track, layerOrient, adjLayer)]
                if pts:
                    via_points[net_name][pin_name][layer] = pts

    obsts = {}
    markUnusedPins(nets, insts, pins, obsts)

    net_rects = defaultdict(list)
    for net in nets:
        for layer, rects in net._sol.items():
            net_rects[layer].extend(rects)

    rtree_cache = buildTree(nets, insts, obsts)
    net_vertex_cache = {net._name: precompute_net_vertices(net._name, guide_data, track, layerOrient, rtree_cache) for net in nets}

    return track, net_subnets, via_points, rtree_cache, net_vertex_cache, netDict

def route_nets(nets, net_subnets, via_points, net_vertex_cache, rtree_cache, insts, obsts):
    #from collections import defaultdict

    all_net_paths = {}
    

    for i, net in enumerate(nets):
        print(f"Routing net {i+1}/{len(nets)}: {net._name}")
        all_paths = []
        subnets = net_subnets[net._name]
        if not subnets:
            continue
        vertices_set = net_vertex_cache[net._name]
        blocked = filter_vertices_by_rtree(vertices_set, rtree_cache, net._name)
        routed_vertices_so_far = set()

        for pin1, pin2 in subnets:
            pin1_vias = via_points[net._name].get(pin1, {})
            pin2_vias = via_points[net._name].get(pin2, {})
            sources = [(x, y, l) for l, pts in pin1_vias.items() for (x, y, l) in pts] + list(routed_vertices_so_far)
            targets = [(x, y, l) for l, pts in pin2_vias.items() for (x, y, l) in pts]
            if not sources or not targets:
                continue

            unblocked_vertices = [v for v in vertices_set if v not in blocked]
            G = build_3d_grid_graph(unblocked_vertices, layerOrient, via_points[net._name])
            path = try_astar_with_virtuals(G, sources, targets)

            if path is None:
                path = force_route(net._name, vertices_set, blocked, sources, targets, track, layerOrient, adjLayer, guide_data, via_points)
            if path:
                print("routed")
                all_paths.append(path[::-1])
                expand_sources_with_path(path[1:-1], routed_vertices_so_far)

        all_net_paths[net._name] = all_paths
        for path in all_paths:
            rects = convert_path_to_rects(path[1:-1], layerWidth)
            for layer, rs in rects.items():
                if layer not in net._sol:
                    net._sol[layer] = []
                net._sol[layer].extend(rs)

        rtree_cache = buildTree(nets, insts, obsts)

    return all_net_paths

def write_solution(nets, ideff, output_path, all_net_paths):
    write_paths_to_def(
        d=ideff,
        output_path=output_path,
        all_net_paths=all_net_paths,
        layer_width=layerWidth
    )
    
def reorder_nets_by_priority(all_nets, high_priority_names):
    high_priority = [net for net in all_nets if net._name in high_priority_names]
    low_priority = [net for net in all_nets if net._name not in high_priority_names]
    return  high_priority + low_priority  
    
def remove_rects_from_rtree(rtree_cache, net_id, net_sol):
    """
    Removes all entries from rtree_cache that match the given net_id and overlap with net_sol.
    """
    for layer, rects in net_sol.items():
        if layer not in rtree_cache:
            continue
        idx = rtree_cache[layer]
        for r in rects:
            box = (r.ll.x, r.ll.y, r.ur.x, r.ur.y)
            matches = list(idx.intersection(box, objects=True))
            for match in matches:
                if match.object == net_id:
                    try:
                        idx.delete(match.id, box)
                    except Exception as e:
                        print(f"Warning: failed to remove rect {box} from layer {layer}: {e}")
           
           
if __name__ == '__main__':
    # [Argument parsing remains the same]
    import argparse
    import time
    ap = argparse.ArgumentParser(description="Detailed Router Phases 1-7")
    ap.add_argument('-l','--leff', required=True)
    ap.add_argument('-i','--ideff', required=True)
    ap.add_argument('-g','--guide', required=True)
    ap.add_argument('-o','--output', required=True)
    ap.add_argument("-p", "--plot", action='store_true')
    args = ap.parse_args()
    
 
    start = time.time()

    leff, ideff, guide_data, insts, pins, nets = load_design_files(args.leff, args.ideff, args.guide)
    track, net_subnets, via_points, rtree_cache, net_vertex_cache, netDict = preprocess_routing_data(ideff, nets, insts, pins, guide_data)
    all_net_paths = route_nets(nets, net_subnets, via_points, net_vertex_cache, rtree_cache, insts, obsts={})
    write_solution(nets, ideff, args.output, all_net_paths)
    end = time.time()
    print(f"\n Routing completed in {end - start:.2f} seconds")
    # Post-routing analysis
    netDict, insts, pins, obsts, nets, spacingnets, openset = loadAndCheck(args.output, args.ideff, args.leff, args.plot)
    violating_nets = spacing_net(spacingnets)
    #print(f"Violating nets: {violating_nets}")
    #print(f"Open nets: {openset}")

    high_priority_names = violating_nets | openset
    nets = reorder_nets_by_priority(nets, high_priority_names)
    # print(f"Prioritized net order: {[net._name for net in nets]}")
    # # for net in nets:
    # #     if net._name in high_priority_names:
    # #         remove_rects_from_rtree(rtree_cache, net._id,net._sol)
    # #         net._sol.clear()

    start_rrr = time.time()
    #RRR    
    rrr(iter = 2,
    open_nets=openset,
    violating_nets=violating_nets,
    netDict=netDict,
    net_subnets=net_subnets,
    all_net_paths=all_net_paths,
    insts=insts,    
    obsts=obsts,
    nets=nets,    
    guide_data=guide_data,    
    layerOrient=layerOrient,
    via_points=via_points,
    layerwidth=layerWidth,
    
    output_def=args.output,
    original_def=args.ideff,
    lef=args.leff
)
    end_time = time.time()    
    print(f"total_time={end_time-start}")
    
    