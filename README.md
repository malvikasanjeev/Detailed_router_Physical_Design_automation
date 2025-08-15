# Detailed Router
# VLSI Detailed Router in Python

This project is a Python-based detailed router developed for the EE5333: Physical Design Automation course. Its primary focus is to achieve **complete net connectivity** by generating detailed routes that result in **zero open nets**. The core of the router uses an **A\* search algorithm** for pathfinding and a powerful **Rip-up and Reroute (RRR)** strategy to iteratively resolve connectivity issues.



## 📝 Description

In VLSI physical design, a detailed router's job is to find the exact metal paths for every connection. While minimizing Design Rule Checking (DRC) violations is important, the most critical goal is ensuring every net is fully connected.

This program automates the routing process with a focus on connectivity by:
1.  Parsing standard industry files: **LEF** (Library Exchange Format), **DEF** (Design Exchange Format), and a **GUIDE** file for global routing hints.
2.  Implementing a cost-based A\* search to find an initial path for each connection.
3.  Employing an aggressive Rip-up and Reroute (RRR) mechanism that prioritizes and fixes any nets that fail to connect, iterating until there are zero open nets.
4.  Generating a final `output.def` file containing the complete, fully connected routing solution.

## ✨ Features

* **Connectivity-Focused Routing:** The entire algorithm is optimized to ensure every part of every net is successfully connected.
* **A\* Search Algorithm:** Finds an efficient initial path on the routing grid, navigating around obstacles and existing routes.
* **Rip-up and Reroute (RRR):** The key feature for achieving 100% connectivity. It iteratively identifies and re-routes failed nets to resolve open connections.
* **LEF/DEF/GUIDE Parser:** Reads and interprets standard VLSI design files to understand the routing problem.
* **Validation Checker (`checker.py`):** A standalone utility to analyze a routed DEF file, confirming the absence of open nets and reporting any remaining spacing violations.
* **Visualization:** Includes options to plot the final routing using `matplotlib` for debugging and analysis.

## 📂 Project Structure

The project is organized into the following directories:

```
.
├── def/              # Input DEF files for different designs
├── guide/            # Input GUIDE files for global routing
├── lef/              # Input LEF file for library and layer info
├── output/           # Output directory for routed DEF files
└── src/              # All Python source code
    ├── detailed_router.py  # The main routing script
    ├── checker.py          # The validation and checking script
    
```



## 🧠 Algorithm 

This router employs a multi-stage process to achieve its primary goal of 100% net connectivity. The core logic is built around the A\* search algorithm, followed by a powerful iterative repair phase known as Rip-up and Reroute (RRR).

### 1. Initialization and Grid Formation
* **Parsing:** The router begins by parsing the input **LEF** and **DEF** files to build an internal model of the chip's physical layout. This includes defining the properties of each metal layer (direction, spacing, width) and identifying all physical obstacles (macros, pre-existing routes).
* **Cost Map:** A 3D cost grid representing the chip's layers is generated. Each cell in the grid is assigned a cost, with obstacles being marked as infinitely expensive (unroutable).

### 2. A\* Pathfinding
* **Core Routing Engine:** For each two-point connection (a "pin-pair") within a net, the router uses an **A\* search algorithm** to find the optimal path.
* **Heuristic Function:** The A\* algorithm is guided by a heuristic that estimates the cost to the destination. This ensures it intelligently explores promising paths first, making it far more efficient than a simple brute-force search.
* **Cost Management:** The cost function considers:
    * **Path Length:** The physical distance of the route.
    * **Via Cost:** A penalty for changing metal layers, which is crucial for creating efficient routes.
    * **Congestion:** As routes are added, the costs of the grid cells they occupy are increased. This discourages future routes from using the same congested areas, naturally preventing shorts.

### 3. Rip-up and Reroute (RRR)
After the initial routing pass, the router enters the critical RRR phase to resolve any remaining open (unconnected) nets.
* **Validation:** The router first identifies all nets that failed to connect completely.
* **Iterative Repair:** It then enters a loop where it systematically:
    1.  **Rips Up:** Removes the failed routes from the grid, resetting the costs of the cells they previously occupied.
    2.  **Reroutes:** Attempts to route these failed nets again. Because other nets have already been routed, the cost landscape of the chip is different, often opening up new potential solutions.
* **Forced Routing:** If a net continues to fail even during RRR, the router employs aggressive "forced" techniques. This includes expanding the valid routing area (the "guide region") or elevating parts of the route to different layers to bypass tricky obstacles, ensuring a connection is made.

This systematic, multi-stage approach guarantees robust performance, allowing the router to successfully find a fully connected solution even in complex and congested designs.

## 🏆 Performance & Results

The router was benchmarked against several standard circuit designs to validate its primary objective: achieving **zero open nets**. The final model incorporates an iterative Rip-up and Reroute (RRR) strategy with "forced routing" techniques, such as expanding the guide region and unblocking paths, to ensure full connectivity.

The table below summarizes the performance on the final, most successful configuration.

| Model | # Nets | Open Nets | DRCs | Time (s) |
| :--- | :--- | :--- | :--- | :--- |
| c17 | 22 | **0** | 5 | 1.7 |
| add5 | 60 | **0** | 56 | 7.5 |
| c432 | 197 | **0** | 325 | 64 |
| spm | 307 | **0** | 707 | 508 |
| c499 | 362 | **0** | 645 | 226 |
| c6288 | 1525 | **0** | 2879 | 2917 |

### Key Takeaways:

* [cite_start]**100% Connectivity Achieved:** The final model successfully routed all nets for every benchmark circuit, resulting in **zero open nets**[cite: 34]. This was the primary goal of the project.
* [cite_start]**Scalability:** The router demonstrates it can handle complex designs, scaling from the simple `c17` (22 nets) to the much larger `c6288` (1525 nets)[cite: 34].
* [cite_start]**Trade-off Management:** While the focus was on connectivity, the router still produces a solution with a quantifiable number of DRC violations[cite: 34]. In a real-world scenario, this output would be the input for a further DRC cleanup stage.

[cite_start]The initial models without this "forced routing" strategy often left several nets open[cite: 30]. [cite_start]For example, the base model left `c432` with 7 open nets and `c6288` with 33 open nets[cite: 30]. This comparison highlights the critical success of the final RRR algorithm.