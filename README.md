# GRASP: Graph Read Aligner with Scalable Partitioning

**GRASP** is a novel, high-performance tool designed for aligning sequencing reads to massive pangenome graphs. Built as a significant evolution over the existing **PanAligner**, GRASP introduces a scalable parallelization framework that overcomes the memory and performance bottlenecks of monolithic aligners.

By leveraging hierarchical graph partitioning and distributed processing, GRASP enables the analysis of complex, terabyte-scale pangenomes on standard hardware, making pangenomics accessible to a wider research community.

## 🚀 Key Features

*   **Scalable Partitioning:** Uses **METIS** to intelligently divide large graphs into balanced, manageable subgraphs (partitions), drastically reducing memory footprint.
*   **Parallel Execution:** Distributes the alignment workload across multiple cores or nodes. Each partition is processed independently.
*   **Halo Region Support:** Automatically identifies and duplicates "halo" (boundary) nodes to preserve alignment context across partition cuts, ensuring high accuracy.
*   **Smart Read Dispatching:** A minimizer-based global index routes reads to the specific partition(s) most likely to contain the correct alignment.
*   **Standard Output:** Produces alignments in the standard **GAF (Graph Alignment Format)**.

## 🛠️ Prerequisites

Before building GRASP, ensure you have the following dependencies installed:

*   **C++ Compiler:** GCC (g++) supporting C++20 standard.
*   **Zlib:** Compression library (`zlib1g-dev`).
*   **METIS:** Serial Graph Partitioning and Fill-reducing Matrix Ordering (`libmetis-dev`).
*   **GKlib:** A library used by METIS (`libgklib`).
*   **GNU Make:** For building the project.

## 📦 Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/KshitijThareja/GRASP.git
    cd GRASP
    ```

2.  **Build the project:**
    ```bash
    make
    ```
    This will compile the executable (currently named `PanAligner`).

## 💻 Usage

GRASP supports two modes of operation: **Standard (Monolithic)** and **Parallel (Partitioned)**.

### 1. Parallel Mode (Recommended for Large Graphs)
This mode automatically partitions the graph, dispatches reads, runs parallel alignments, and merges the results.

**Using the Automation Script:**
The easiest way to run in parallel mode is using the provided script:

```bash
./run_parallel.sh <GRAPH_GFA> <READS_FASTQ> [THREADS]
```

*   `<GRAPH_GFA>`: Input pangenome graph in GFA format (e.g., `graph.gfa`).
*   `<READS_FASTQ>`: Input reads in FASTA or FASTQ format.
*   `[THREADS]`: Number of parallel threads/jobs to launch (default: 4).

**Example:**
```bash
./run_parallel.sh yeast_pangenome.gfa reads.fq 8
```

**Output:**
*   `final_alignment.gaf`: The merged alignment results.
*   `partition_*.gfa`: Generated partition subgraphs (intermediate).
*   `reads_part_*.fq`: Dispatched reads per partition (intermediate).

---

### 2. Standard Mode
For smaller graphs that fit in memory, you can run GRASP as a standard single-process aligner.

```bash
./PanAligner [options] <target.gfa> <query.fa|fq> > output.gaf
```

**Options:**
*   `-t INT`: Number of threads (default: 4).
*   `-c`: Generate CIGAR string in output.
*   `--partition`: (Internal flag) Triggers the partitioning and dispatching phase.

## 🧩 Architecture Overview

GRASP operates in three stages:

1.  **Partitioning & Dispatching:**
    *   The global graph is analyzed for Weakly Connected Components (WCCs).
    *   Large components are split using the METIS k-way partitioning algorithm.
    *   A global minimizer index is built to map sequence features to partition IDs.
    *   Incoming reads are "voted" into the appropriate partition(s).

2.  **Parallel Alignment:**
    *   Independent worker processes load specific partitions (`partition_X.gfa`) and their assigned reads (`reads_part_X.fq`).
    *   Alignment is performed using a seed-chain-extend approach, optimized for the local subgraph.

3.  **Merging:**
    *   Results from all workers are aggregated into a single, consistent GAF file.

## 📄 Output Format
The output follows the **GAF (Graph Alignment Format)** specification:
```text
QueryName  QueryLen  QueryStart  QueryEnd  Strand  Path  PathLen  PathStart  PathEnd  Matches  BlockLen  MapQ  ...
```
