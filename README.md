# Temporal Stable Community Search (TSCS)

## Introduction

This project addresses the problem of temporal stable community search (TSCS). It includes both an online TSC search algorithm and an index-based TSC search algorithm, providing efficient solutions to find stable communities on temporal graphs.

## Environment

The TSCS codes are implemented and tested in the following development environment:

- **Hardware**: Intel(R) Xeon(R) Gold 2.40GHz CPU and 384GB of memory
- **Operating System**: Ubuntu 20.04.4 LTS (GNU/Linux 5.13.0-40-generic x86_64)
- **C++ Version**: 17

## Datasets

The dataset files can be found on [Network Repository](https://networkrepository.com/network-data.php) and [Konect](http://konect.cc/networks/). Detailed statistics of these datasets are provided in our paper.

## How to Run the Codes

### Online TSC Search

The program requires three arguments: the path to the graph file, the path to the query file, and the path to the answer file.

#### Code Compilation

```bash
chmod +x build_OTSCS.sh
./build_OTSCS.sh
```

#### Run code

```
./OTSCS [Input_graph_path] [Input_query_path] [Output_answer_path]
```

### Index-based TSC search

The program requires at least three arguments: the path to the graph file, the path to the query file, the path to the answer file, the option on AG-index or ASF-index, and the option on AITE or BITE.

Additional functions like `load index` and `save index` are supported. Note that AG-index and ASF-index share a part of the index, referred to as the core index. The ASF-index and AG-index are divided into two parts: the corresponding core index and the rest part, to reduce redundant calculations.

#### Code compilation

```bash
chmod +x build_TSCS.sh
./build_TSCS
```

#### Run code

##### Basic usage

```bash
./TSCS [Input_graph_path] [Input_query_path] [Output_answer_path] -AG/-ASF -AITE/-BITE 
```

##### Options

- `-ASF` or `-AG`: Use ASF-index or AG-index.
- `-AITE` or `-BITE`: Use AITE or BITE for query. 
- `-sc <path>`: Save the core index (the shared part of ASF/AG-index) to the specified path.
- `-lc <path>`: Load the core index (the shared part of ASF/AG-index) from the specified path.
- `-st <path>`: Save the ASF/AG index (the rest part, core index not included) to the specified path.
- `-lt <path>`: Load the ASF/AG index (rest part, core index not included) from the specified path.

##### Example

```
./TSCS graph.txt queries.txt answers.txt -ASF -BITE -sc core_index.txt -st asf_index.txt
```

This command will:

1. Use ASF-index and BITE as the query algorithm.
2. Save the core index to `core_index.txt` (`-sc core_index.txt`).
3. Save the rest part of ASF index to `asf_index.txt` (`-st asf_index.txt`).


```
./TSCS graph.txt queries.txt answers.txt -AG -AITE -lc core_index.txt -st ag_index.txt
```

This command will:

1. Use ASF-index and AITE as the query algorithm.
2. Load the core index from `core_index.txt` (`-lc core_index.txt`).
3. Save the rest part of AG index to `ag_index.txt` (`-st ag_index.txt`).

### Demo
A demo graph and the corresponding query can be found as `demo_graph.txt` and `demo_query.txt`.

- **`demo_graph.txt`**: 
  - The first line contains three integers: $n$, $m$, and $t_{max}$ respectively.
  - The following $m$ lines represent the $m$ edges in the graph. Each line contains three integers $u$, $v$, and $t_s$ representing an undirected edge between nodes $u$ and $v$ at time $t_s$.

- **`demo_query.txt`**:
  - The first line contains a single integer representing the number of queries.
  - Each subsequent line contains four integers: $t_s$, $t_e$, $q$, and $k$, representing a TSC search.


After compiling TSCS, you can test it using the following command:

```bash
./TSCS demo_graph.txt demo_query.txt answers.txt -AG -AITE
```

This command will solve the TSC search problem using the AG-index and AITE algorithm on the demo graph and queries, outputting the result to `answers.txt`.