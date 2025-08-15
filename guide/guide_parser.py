class GuideFileParser:
    def __init__(self, file_path):
        """Initializes the parser with the guide file path"""
        self.file_path = file_path
        self.guide_data = {}

    def parse(self):
        """Parses the guide file into a dictionary of nets and segments"""
        current_net = None

        with open(self.file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                
                parts = line.split()
                
                # If line has 5 parts and first 4 are digits, it's a segment
                if len(parts) == 5 and all(part.isdigit() for part in parts[:4]):
                    if current_net is None:
                        print(f"Warning: Segment found without a net name: {line}")
                        continue
                    
                    x1, y1, x2, y2, layer = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), parts[4]
                    self.guide_data[current_net].append((x1, y1, x2, y2, layer))
                else:  # It's a net name
                    current_net = line
                    if current_net not in self.guide_data:
                        self.guide_data[current_net] = []

        return self.guide_data

    def get_net_name(self, index=0):
        """Returns a net name as a string at the given index"""
        net_names = list(self.guide_data.keys())
        if index < len(net_names):
            return net_names[index]
        return None

    def get_segment_coordinates(self, net_name, segment_index=0):
        """
        Returns the coordinates for a specific net segment
        
        Returns:
            tuple: (x_ll, y_ll, x_ur, y_ur) for the specified segment
        """
        segments = self.guide_data.get(net_name, [])
        if segments and segment_index < len(segments):
            return segments[segment_index][:4]  # First 4 elements are coordinates
        return None
        
    def get_segment_layer(self, net_name, segment_index=0):
        """
        Returns layer name as string for a specific segment
        """
        segments = self.guide_data.get(net_name, [])
        if segments and segment_index < len(segments):
            return segments[segment_index][4]  # Layer is the 5th element
        return None
        
    def get_segment_details(self, net_name, segment_index=0):
        """
        Returns all details for a specific segment
        
        Returns:
            tuple: (net_name, x_ll, y_ll, x_ur, y_ur, layer)
        """
        segments = self.guide_data.get(net_name, [])
        if segments and segment_index < len(segments):
            x_ll, y_ll, x_ur, y_ur, layer = segments[segment_index]
            return (net_name, x_ll, y_ll, x_ur, y_ur, layer)
        return None

    def get_net_names(self):
        """Returns list of all net names"""
        return list(self.guide_data.keys())

    def get_segments(self, net_name):
        """Returns all segments for a given net"""
        return self.guide_data.get(net_name, [])

    def get_layers(self):
        """Returns set of all unique layer names"""
        layers = set()
        for segments in self.guide_data.values():
            for segment in segments:
                layers.add(segment[4])
        return layers

    def get_segment_count(self, net_name):
        """Returns the number of segments for a given net"""
        return len(self.guide_data.get(net_name, []))
    
    
if __name__ == "__main__":
    guide_file_path = "/home/malavika1606/EE5333_projects/guide/c17.guide"
    parser = GuideFileParser(guide_file_path)
    parser.parse()
    
    # Get a specific segment's details
    net = "N1"
    segment_details = parser.get_segment_details(net)
    if segment_details is None:
        print(f"No segments found for net inside main: {net}")
    else:
        net_name, x_ll, y_ll, x_ur, y_ur, layer = segment_details
        print(f"Net: {net_name}, Coords: ({x_ll},{y_ll}) to ({x_ur},{y_ur}), Layer: {layer}")