# VLSI Detailed Router

A Python-based detailed router for VLSI physical design that guarantees **100% net connectivity** using an A* search algorithm and a Rip-up and Reroute (RRR) strategy. Developed for the EE5333: Physical Design Automation course.

## Performance & Results

The router was benchmarked against several standard circuit designs to validate its primary objective of achieving **zero open nets**.

| Model | # Nets | Open Nets | DRCs | Time (s) |
| :--- | :--- | :--- | :--- | :--- |
| c17 | 22 | **0** | 5 | 1.7 |
| add5 | 60 | **0** | 56 | 7.5 |
| c432 | 197 | **0** | 325 | 64 |
| spm | 307 | **0** | 707 | 508 |
| c499 | 362 | **0** | 645 | 226 |
| c6288 | 1525 | **0** | 2879 | 2917 |

#### Key Takeaways:
* **Guaranteed Connectivity:** Successfully routed all nets for every benchmark, resulting in **zero open nets**.
* **Scalability:** Proven to handle complex designs, scaling from 22 nets (`c17`) to 1525 nets (`c6288`).
* **Trade-off Management:** While prioritizing connectivity, the router produces a solution with a quantifiable number of DRC violations, suitable for a subsequent cleanup stage.

##  Key Features

* **A\* Search Algorithm:** Finds efficient initial paths on the routing grid while navigating obstacles.
* **Rip-up and Reroute (RRR):** An aggressive, iterative repair mechanism that guarantees 100% connectivity by fixing open nets.
* **Forced Routing:** Employs advanced techniques like expanding the search area to ensure difficult connections are completed.
* **Industry Standard I/O:** Parses LEF/DEF files for input and generates a fully routed DEF file as output.
* **Validation Checker:** Includes a `checker.py` utility to confirm connectivity and report on DRC violations.

## How It Works

The router follows a systematic, multi-stage process:

1.  **Initialization:** Parses LEF/DEF files to build a 3D cost grid of the chip's layers and obstacles.
2.  **A\* Pathfinding:** Routes each connection using the A\* algorithm, which finds the lowest-cost path considering distance, vias, and congestion.
3.  **Iterative Repair (RRR):** After the first pass, it identifies any open nets, rips up their failed routes, and intelligently reroutes them until connectivity is 100% complete.

The project source code is located in the `/src` directory, with benchmark files in `/lef`, `/def`, and `/guide`.

## Technologies Used

* **Language:** Python
* **Libraries:** NumPy, Matplotlib, Rtree
