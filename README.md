# GRASP: Graph Read Aligner with Scalable Partitioning

**GRASP** is a novel, high-performance tool designed for aligning sequencing reads to massive pangenome graphs. Built as a significant evolution over the existing **PanAligner**, GRASP introduces a scalable parallelization framework that overcomes the memory and performance bottlenecks of monolithic aligners.

By leveraging hierarchical graph partitioning and distributed processing, GRASP enables the analysis of complex, terabyte-scale pangenomes on standard hardware, making pangenomics accessible to a wider research community.

## Key Features

*   **Scalable Partitioning:** Uses **METIS** to intelligently divide large graphs into balanced, manageable subgraphs (partitions), drastically reducing memory footprint.
*   **Parallel Execution:** Distributes the alignment workload across multiple cores or nodes. Each partition is processed independently.
*   **Halo Region Support:** Automatically identifies and duplicates "halo" (boundary) nodes to preserve alignment context across partition cuts, ensuring high accuracy.
*   **Smart Read Dispatching:** A minimizer-based global index routes reads to the specific partition(s) most likely to contain the correct alignment.
*   **Standard Output:** Produces alignments in the standard **GAF (Graph Alignment Format)**.

## Prerequisites

Before building GRASP, ensure you have the following dependencies installed:

*   **C++ Compiler:** GCC (g++) supporting C++20 standard.
*   **Zlib:** Compression library (`zlib1g-dev`).
*   **METIS:** Serial Graph Partitioning and Fill-reducing Matrix Ordering (`libmetis-dev`).
*   **GKlib:** A library used by METIS (`libgklib`).
*   **GNU Make:** For building the project.

## Installation

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

## Usage

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
*   `-x lr`: Long-read preset (recommended for reads ≥1 kbp; sets appropriate minimizer and chaining parameters).
*   `--partition`: Triggers the graph partitioning and read dispatching phase.
*   `--nparts INT`: Number of partitions K for METIS graph partitioning (default: 4). Higher K reduces per-worker memory but increases partition overhead.
*   `--halo-depth INT`: BFS halo expansion depth d around partition boundaries (default: 1). Halo nodes are duplicated into neighbouring partitions to preserve alignment context at cuts.

## Architecture Overview

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

## Output Format
The output follows the **GAF (Graph Alignment Format)** specification:
```text
QueryName  QueryLen  QueryStart  QueryEnd  Strand  Path  PathLen  PathStart  PathEnd  Matches  BlockLen  MapQ  ...
```

---

## Benchmarking & Evaluation

GRASP ships with three tools for reproducing the paper's experiments.

### `simulate_reads.py`
Generates simulated long reads from a reference FASTA. Read names encode the true origin (`<seqname>_<start>_<end>`) so alignment accuracy can be verified later.

```bash
python3 simulate_reads.py \
    --fasta  reference.fna \
    --output reads.fasta   \
    --length 10000         \
    --num    3000          \
    --seed   42
```

| Option | Description |
|--------|-------------|
| `--fasta` | Input reference FASTA |
| `--output` | Output FASTA path |
| `--length INT` | Read length in bp (default: 10000) |
| `--num INT` | Number of reads to generate (default: 5000) |
| `--step INT` | Step between tiled reads; 0 = random sampling (default: 0) |
| `--error FLOAT` | Per-base substitution error rate 0–1 (default: 0) |
| `--seed INT` | Random seed (default: 42) |

---

### `eval_accuracy.py`
Evaluates alignment accuracy in two modes.

**Agreement mode** — compares GRASP output against a PanAligner baseline; reports the fraction of reads whose primary alignment path matches:
```bash
python3 eval_accuracy.py \
    --mode     agreement \
    --baseline baseline_panaligner.gaf \
    --query    grasp_output.gaf
```

**Position mode** — for simulated reads, checks whether each aligned path covers the read's true origin using SN/SO tags from the GFA:
```bash
python3 eval_accuracy.py \
    --mode  position \
    --gfa   pangenome.gfa \
    --query grasp_output.gaf
```

---

### `benchmark.sh`
Orchestrates the full set of paper experiments across two datasets (yeast and Drosophila).

```bash
bash benchmark.sh [yeast|drosophila|all] [fig1|fig2|fig3|fig4|all]
```

| Figure | What it measures |
|--------|-----------------|
| Fig 1 | Alignment accuracy and cut-edge fraction vs K (# partitions) |
| Fig 2 | Alignment accuracy vs halo depth d |
| Fig 3 | Peak RAM and wall time: PanAligner vs GRASP Full vs GRASP Align-Only |
| Fig 4 | Scalability (speedup, RAM, time) vs K on Drosophila |

Results are written as TSV files to `results/`.

Key environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `FIG3_K` | `4` | K used for the Fig 3 resource comparison run |
| `FIG4_K_VALUES` | `"1 2 4 6"` | K values swept in the Fig 4 scalability experiment |

---

### `run_sequential.sh`
Memory-safe wrapper around `benchmark.sh`. Runs each case one at a time to avoid out-of-memory crashes on large graphs. Accuracy figures (Fig 1, Fig 2) are run with serialised alignment workers; resource figures (Fig 3, Fig 4) keep parallel workers so timing measurements remain valid.

```bash
bash run_sequential.sh [yeast|drosophila|all] [fig1|fig2|fig3|fig4|all]
```

Override K limits via the same environment variables as `benchmark.sh`:
```bash
FIG3_K=4 FIG4_K_VALUES="1 2 4 6" bash run_sequential.sh all all
```
