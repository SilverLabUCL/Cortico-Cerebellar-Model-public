# [PAPER TITLE]: simulation and analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

C++ simulation code and MATLAB analysis scripts accompanying the manuscript
**"[PAPER TITLE]"** by [AUTHORS] ([preprint link / DOI], submitted to [JOURNAL]).

> **Maintainer notes (delete before publishing).** Replace every `[PLACEHOLDER]` and `TODO`.
> Items marked **[CONFIRM]** are my best reading of your notes.
>
> <details><summary>Open items checklist</summary>
>
> - [ ] Overview paragraph and model descriptions (Section 1)
> - [ ] List all `--learnmode` values and write parameter descriptions (Section 5, `docs/PARAMETERS.md`)
> - [ ] Plain-text output file names and the location of each JSON file (Section 7)
> - [ ] Meaning of the `'parameters'` option and full option lists of the MATLAB plotting scripts (Section 8)
> - [ ] How `read_simulations_release.m` distinguishes single simulations from sweeps (Section 8)
> - [ ] Structure of `figure_sim_commands.txt` and the figure table (Section 9)
> - [ ] Zenodo DOI, contents, size, folder layout and where to unpack (Sections 9 and 10)
> - [ ] Peak memory use (Section 11)
> - [ ] Dockerfile and Docker instructions (Section 13)
> - [ ] Paper title, authors, DOIs, BibTeX, contact, funding (Sections 1, 16, 17)
> - [ ] Verify the from-source install steps for GSL and Armadillo (Section 3)
> - [ ] Add `CITATION.cff` and tag the release (e.g. `v1.0-paper`)
> - [ ] Optional / low priority: reference output and check script (Section 6)
> </details>

---

## Contents

1. [Overview](#1-overview)
2. [Repository structure](#2-repository-structure)
3. [Requirements](#3-requirements)
4. [Compilation](#4-compilation)
5. [Running the program](#5-running-the-program)
6. [Test simulations](#6-test-simulations)
7. [Output files](#7-output-files)
8. [MATLAB workflow](#8-matlab-workflow)
9. [Reproducing the paper](#9-reproducing-the-paper)
10. [Data (Zenodo)](#10-data-zenodo)
11. [Computational resources and runtime](#11-computational-resources-and-runtime)
12. [Reproducibility notes](#12-reproducibility-notes)
13. [Docker](#13-docker)
14. [Troubleshooting](#14-troubleshooting)
15. [License](#15-license)
16. [How to cite](#16-how-to-cite)
17. [Contact and acknowledgements](#17-contact-and-acknowledgements)

---

## 1. Overview

[3-5 sentences: the scientific question the paper addresses, and what the code does.
Do not repeat the methods; refer to the paper.]

The repository contains:

- two C++ programs that simulate [DESCRIBE MODEL CLASSES]:
  - `RNN_only.cpp`: [one-line description of the architecture]
  - `RNN+FFL.cpp`: [one-line description; spell out FFL]
- MATLAB scripts that load the simulation output, compute statistics and generate the figures,
- a text file (`figure_sim_commands.txt`) listing the exact commands used for the simulations in the paper,
- a Dockerfile for a reproducible build environment.

The repository does **not** contain the raw simulation output of the paper. The processed data
behind each figure are archived on Zenodo (see [Data](#10-data-zenodo)).

## 2. Repository structure

The C++ programs and the MATLAB scripts live in the same (top-level) folder.

```
REP/
├── RNN_only.cpp                  # simulation: RNN only
├── RNN+FFL.cpp                   # simulation: RNN + FFL
├── headers/                      # C++ header files (including json.hpp from nlohmann/json)
├── simulation_results/           # all output files (written by the programs)
├── read_simulations_release.m    # MATLAB: load simulation output into a struct
├── fig_activity_release.m        # MATLAB: plots for single simulations
├── fig_RNN_scan_release.m        # MATLAB: analysis + plots, RNN_only parameter sweeps
├── fig_RNN_FFL_scan_release.m    # MATLAB: analysis + plots, RNN+FFL parameter sweeps
├── fig_multi_scan_release.m      # MATLAB: analysis + plots combining several structs
├── findIntersect.m               # MATLAB: helper (stability margins from Floquet analysis)
├── figure_sim_commands.txt       # commands used for the simulations in the paper
├── Dockerfile
├── LICENSE
└── README.md
```

`simulation_results/` contains a `.gitkeep` file so that the (otherwise empty) folder
exists after cloning.

## 3. Requirements

### Tested environment

| Component | Version |
|---|---|
| Operating system | Ubuntu 24.04 (Windows is not supported) |
| Compiler | g++ 13.3.0 |
| [Armadillo](https://arma.sourceforge.net/) | 14.4.0 (built from source) |
| [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) | 2.8 (built from source) |
| [Boost](https://www.boost.org/) | 1.83.0 (components used: `program_options`, `version`) |
| BLAS / LAPACK | pthread-enabled OpenBLAS (used for both, accessed through Armadillo) |
| OpenMP | provided by g++ (`-fopenmp`) |
| [nlohmann/json](https://github.com/nlohmann/json) | 3.11.2 (bundled as `headers/json.hpp`, nothing to install) |
| MATLAB | R2024b (24.2.0.3212159), Update 9 |
| MATLAB toolboxes | Statistics and Machine Learning Toolbox |

MATLAB scripts read the JSON output with `jsondecode`.

### Installing the dependencies

**Compiler, Boost and OpenBLAS** (Ubuntu 24.04 packages; the packaged Boost is 1.83.0):

```bash
sudo apt update
sudo apt install build-essential libboost-program-options-dev libboost-dev libopenblas-dev liblapack-dev cmake
```

**GSL 2.8 and Armadillo 14.4.0** were built from source. Download the source tarballs from the
[GSL](https://www.gnu.org/software/gsl/) and [Armadillo](https://arma.sourceforge.net/download.html)
websites, then, for example:

```bash
# GSL
tar xzf gsl-2.8.tar.gz && cd gsl-2.8
./configure && make -j && sudo make install
cd ..

# Armadillo
tar xf armadillo-14.4.0.tar.xz && cd armadillo-14.4.0
cmake . && make -j && sudo make install
cd ..

sudo ldconfig
```

> **TODO [CONFIRM]:** these are the standard build steps; adjust them to match how you actually
> built the libraries (install prefix, configure options). If a different install prefix is used,
> add the corresponding `-I`/`-L` flags to the compile commands. The
> [Docker image](#13-docker) pins the exact versions.

## 4. Compilation

Compile from the repository root:

```bash
g++ -O3 -o RNN_only.out RNN_only.cpp -lgsl -lgslcblas -lm -larmadillo -lboost_program_options -fopenmp

g++ -O3 -o RNN+FFL.out RNN+FFL.cpp -lgsl -lgslcblas -lm -larmadillo -lboost_program_options -fopenmp
```

| Flag | Purpose |
|---|---|
| `-O3` | highest standard optimisation level (important for simulation speed) |
| `-lgsl -lgslcblas -lm` | link GSL, its CBLAS interface and the math library |
| `-larmadillo` | link Armadillo (which uses OpenBLAS/LAPACK underneath) |
| `-lboost_program_options` | link Boost.Program_options (command-line parsing) |
| `-fopenmp` | enable OpenMP multithreading |

Compiling both programs takes about 40 s in total.

Before running either program, set the number of OpenMP threads (the paper simulations used 6):

```bash
export OMP_NUM_THREADS=6
```

The thread count changes only the runtime, not the numerical results.

## 5. Running the program

Both programs take their settings from the command line, parsed with
[Boost.Program_options](https://www.boost.org/doc/libs/release/doc/html/program_options.html).
The options used in the test simulations are:

| Option | Meaning | Values used here |
|---|---|---|
| `--learnmode` | what the program does | `none` (single simulation), `scan_lambda` (parameter sweep over `lambda`) **[CONFIRM: list all modes]** |
| `--Ncomp` | number of frequency components of the target signal | `3` (single simulation), `1` (sweep) |
| `--Nreal` | number of network realisations (repetitions) | `1`, `3`; the paper uses 50 |
| `--seed` | seed of the random number generator (results are deterministic for a given seed) | `-1` |
| `--lambda` | multiplier that scales the frequencies of the target signal | `10` |
| `--Floquet` | compute Floquet multipliers (`ON`/`OFF`) | `OFF` |

The full list of options, their defaults and descriptions is printed by:

```bash
./RNN_only.out --help
./RNN+FFL.out --help
```

> **TODO (later):** parameter descriptions. Suggestion: keep `--help` as the single source of
> truth, keep only the essential options in this table, and put longer descriptions (including
> the contents of `parameters.json` and `flags.json`) in `docs/PARAMETERS.md`, linked from here.

Every run also writes a **log file**, `simulation_results/log_file.dat`, that accumulates a
record of which stage of the program is currently executing (for example `Fourier analysis`).
This is useful for tracking progress in long simulations.
[TODO CONFIRM: is the log file also overwritten by each new run?]

## 6. Test simulations

Run these from the repository root after compiling (and after `export OMP_NUM_THREADS=6`).
All four commands print a short console summary: a few key parameters, the names of the
folders where output is saved, and the Armadillo and Boost versions in use.

| # | Command | Runtime | What it produces |
|---|---|---|---|
| 1 | `./RNN_only.out --learnmode none --Ncomp 3 --Nreal 1 --seed -1 --lambda 10 --Floquet OFF` | ~5 s | A single RNN-only simulation like the one in **Fig. 3b-d** |
| 2 | `./RNN+FFL.out --learnmode none --Ncomp 3 --Nreal 1 --seed -1 --lambda 10 --Floquet OFF` | ~5 s | A single RNN+FFL simulation like the one in **Fig. 3f-h** |
| 3 | `./RNN_only.out --learnmode scan_lambda --Ncomp 1 --Nreal 3 --seed -1 --Floquet OFF` | ~6 min | RNN-only parameter sweep: a low-repetition (3 instead of 50) version of **Fig. 1d** |
| 4 | `./RNN+FFL.out --learnmode scan_lambda --Ncomp 1 --Nreal 3 --seed -1 --Floquet OFF` | ~8 min | RNN+FFL parameter sweep: a low-repetition (3 instead of 50) version of **Fig. 1d** |

Runtimes are for the hardware described in [Section 11](#11-computational-resources-and-runtime).

**Note:** every call of a program **overwrites** the contents of `simulation_results/` from
previous runs. Copy or rename the folder if you want to keep the output of one test before
running the next.

### What the output should look like

All output is written to `simulation_results/` (see [Output files](#7-output-files)).

- **Tests 1 and 2** (single simulations) write the general files (`parameters.json`,
  `flags.json`, population firing rates, `log_file.dat`) and a subfolder `parameter_1/`, which
  also contains `weights_single_trial.json` with the network connectivity.
- **Tests 3 and 4** (sweeps) create, for each parameter set `i` of the sweep, a folder
  `parameter_i/` containing the same kind of files as in the single-simulation case
  **[CONFIRM: whether `weights_single_trial.json` is also written in sweeps]**.

### Checking that a test worked

1. The console summary appears and reports the expected Armadillo (14.4.0) and Boost (1.83.0) versions.
2. The files and folders listed in [Output files](#7-output-files) exist and are non-empty.
3. Load the output in MATLAB with `read_simulations_release.m`, then plot it with
   `fig_activity_release.m` (tests 1 and 2, compare with Fig. 3b-d and 3f-h) or
   `fig_RNN_scan_release.m` / `fig_RNN_FFL_scan_release.m` (tests 3 and 4, compare with Fig. 1d,
   keeping in mind that only 3 repetitions are used). Example calls are given in
   [Section 8](#8-matlab-workflow).

> **Optional (later):** since runs are deterministic for a given seed, a small reference output
> and a check script (comparing within a tolerance) could be added so that users can verify
> their installation automatically.

## 7. Output files

All output is written to `simulation_results/`. Files with the extension `.json` are in JSON
format (written with [nlohmann/json](https://github.com/nlohmann/json)); all other output
files are plain text.

| File / folder | Format | Contents |
|---|---|---|
| `log_file.dat` | text | log of the program's progress |
| `parameters.json` | JSON | all simulation parameters |
| `flags.json` | JSON | all flags |
| [TODO: file names] | text | population firing rates from the test phase and from an open-loop phase |
| `parameter_1/` (`parameter_i/` in sweeps) | | one folder per parameter set, see below |

**Contents of `parameter_1/` (`parameter_i/`)** — for each network realisation (`Nreal`)
**[CONFIRM: assignment of contents to file names]**:

| File | Format | Contents |
|---|---|---|
| `targets.json` | JSON | target signals |
| `stats.json` | JSON | network firing-rate statistics |
| `measures.json` | JSON | network representational measures |
| [TODO: file name] | text | test-phase readout neuron activity |
| `weights_single_trial.json` | JSON | network connectivity (recurrent weights, readout weights, etc.) together with single-unit time constants, gains and thresholds, from a single trial |

Each new run overwrites all content from previous simulations in this folder.

## 8. MATLAB workflow

```
C++ programs  ->  files in simulation_results/  ->  read_simulations_release.m  ->  MATLAB struct
              ->  fig_*_release.m (statistics over realisations + plots)
```

The MATLAB scripts are located in the same folder as the C++ programs. Run them with the
repository root as the current folder. Requires MATLAB R2024b and the Statistics and Machine
Learning Toolbox.

| Script | Purpose |
|---|---|
| `read_simulations_release.m` | Reads the simulation data from the output folder and puts everything into a MATLAB struct. The struct has a slightly different layout for a single simulation and for a parameter sweep. |
| `fig_activity_release.m` | Plots the output of a **single simulation** (e.g. test commands 1 and 2). Can plot principal-component modes and Fourier modes of the population activity. Reads the data from the struct. |
| `fig_RNN_scan_release.m` | Analyses and plots **RNN_only parameter sweeps** (e.g. test command 3). Averages measures extracted from the simulations over multiple realisations and plots, for example, average output errors, principal angles, norms and Floquet multipliers in multiple panels. Reads the struct created by `read_simulations_release.m`. |
| `fig_RNN_FFL_scan_release.m` | Same as `fig_RNN_scan_release.m`, for **RNN+FFL** sweeps (e.g. test command 4). |
| `fig_multi_scan_release.m` | Data analysis and plotting from **several structs**, allowing multiple curves in the same panel. |
| `findIntersect.m` | Helper routine used by the plotting scripts to find the stability margins obtained from the Floquet analysis. |

### Example calls

**Load the simulation output into a struct** (works for a single simulation or a parameter
sweep; the struct is created in the MATLAB workspace):

```matlab
MAT_STRUCT = read_simulations_release('simulation_results/');
```

**Single simulation** (e.g. after test 1 or 2): plot the activity of the readout unit, the
activity of some example neurons in the RNN and the FFL, the principal component modes and,
in this example, the first 6 Fourier modes:

```matlab
fig_activity_release(MAT_STRUCT, 'show output', 'show activity', 'show modes', 6);
```

**RNN_only parameter sweep** (e.g. after test 3): plot statistics of the sweep versus the sweep
parameter. Here: RMSE, task-relevant vector norms and the sine product (as defined in the
paper). `'RMSE log axis'` sets the y-axis of the RMSE panel to log scale and `[1e-3 1]` sets
its y-axis limits:

```matlab
fig_RNN_scan_release(MAT_STRUCT, 'RMSE log axis', [1e-3 1], 'RMSE', 'norms', 'sine product', 'parameters');
```

**RNN+FFL parameter sweep** (e.g. after test 4): same as above, for RNN+FFL data:

```matlab
fig_RNN_FFL_scan_release(MAT_STRUCT, 'RMSE log axis', [1e-3 1], 'RMSE', 'norms', 'sine product', 'parameters');
```

**Several parameter sweeps in the same panels** (here two structs, e.g. one from each program):

```matlab
fig_multi_scan_release('data', MAT_STRUCT_1, MAT_STRUCT_2, 'RMSE log axis', [1e-3 1], 'RMSE', 'norms', 'sine product');
```

> **TODO:** explain the `'parameters'` option, and either list all options of the plotting
> scripts here or point to a header comment in each script. Also state how
> `read_simulations_release.m` distinguishes a single simulation from a sweep
> (automatically, or via an argument).

## 9. Reproducing the paper

The commands used for the simulations in the paper are listed in
**`figure_sim_commands.txt`**, included in this repository.
[TODO: describe how the file is organised, e.g. one block per figure.]

| Figure | Data (Zenodo) | MATLAB script |
|---|---|---|
| Fig. 1 | [file] | [script] |
| Fig. 3 | [file] | [script] |
| ... | ... | ... |
| Fig. 8 | [file] | [script] |

There are two ways to proceed:

- **(a) Regenerate the figures from the processed data (fast).** Download the data from Zenodo
  ([Section 10](#10-data-zenodo)) and run the plotting scripts. No C++ compilation is needed.
- **(b) Rerun the simulations (slow).** Compile the programs, run the relevant commands from
  `figure_sim_commands.txt`, then run the MATLAB scripts on the new output. Expect runtimes of
  up to several hours per figure ([Section 11](#11-computational-resources-and-runtime)).

## 10. Data (Zenodo)

The processed data used to generate the figures are archived at
**[https://doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)** (license: CC BY 4.0).

- **Contents:** processed figure data (text files). Raw simulation output is *not* included.
- **Size:** [TODO]
- **Download and placement:** [TODO: e.g. unpack into `data/` in the repository root, and state
  which path the MATLAB scripts expect]
- The code in this repository is archived on Zenodo through the GitHub integration. Use the
  **version DOI** of the release tagged `[v1.0-paper]` to cite the exact code used in the paper;
  the concept DOI always points to the latest version.

## 11. Computational resources and runtime

All simulations were run on a Lenovo ThinkPad P14s Gen 5 (Intel Core Ultra 7 155H, 22 threads,
32 GB RAM) under Ubuntu 24.04, using 6 OpenMP threads. No cluster was used.

| Task | Approximate runtime |
|---|---|
| Compilation (both programs) | ~40 s |
| Test simulations 1 and 2 | ~5 s each |
| Test simulation 3 | ~6 min |
| Test simulation 4 | ~8 min |
| Typical full parameter sweep (50 repetitions), RNN_only or RNN+FFL, **without** Floquet multipliers | ~1 hour |
| Typical full parameter sweep (50 repetitions), RNN_only or RNN+FFL, **with** Floquet multipliers | ~8 hours |
| Fig. 8 (the most complex architecture, Fig. 8g, including Floquet multipliers) | ~50 hours |

Memory: [TODO: peak RAM usage of a typical run].
Disk space: small; the whole repository including the simulated data is typically smaller than 50 MB.

## 12. Reproducibility notes

- Results are deterministic given a random number seed (`--seed`); `--seed -1` is a fixed seed
  and gives the same output on every run on a given machine.
- The number of OpenMP threads (`OMP_NUM_THREADS`) affects only the runtime, not the numerical results.
- Small numerical differences can occur across machines, compilers or BLAS/LAPACK implementations.
  [TODO: state the expected size of such differences, if known.]

## 13. Docker

A Dockerfile with all dependencies pinned to the versions above is provided.

```bash
docker build -t rep .
docker run --rm -it -v "$(pwd)":/work rep
```

TODO: adapt these commands to the final Dockerfile, and state that MATLAB is not included in the
image (it is proprietary software), so the MATLAB scripts must be run on the host.

## 14. Troubleshooting

- **`cannot find -lgsl` / `-lgslcblas` / `-larmadillo`:** the library is missing or not on the linker path. After a source install into `/usr/local`, run `sudo ldconfig`.
- **`cannot find -lboost_program_options`:** install `libboost-program-options-dev`.
- **Armadillo link or LAPACK errors:** make sure the OpenBLAS/LAPACK development packages are installed.
- **Slow runs (single-threaded):** the program was probably compiled without `-fopenmp`, or `OMP_NUM_THREADS` is not set.
- **Output of an earlier run disappeared:** each run overwrites `simulation_results/`; copy the folder first.
- **Files or folders not found:** run the programs and MATLAB scripts from the repository root.

## 15. License

- **Code:** MIT (see [`LICENSE`](LICENSE)). The bundled `headers/json.hpp` is from
  [nlohmann/json](https://github.com/nlohmann/json) (MIT license).
- **Data on Zenodo:** CC BY 4.0.

## 16. How to cite

If you use this code or data, please cite the paper and the archived code release:

```bibtex
@article{[KEY],
  title   = {[PAPER TITLE]},
  author  = {[AUTHORS]},
  journal = {[JOURNAL / bioRxiv]},
  year    = {[YEAR]},
  doi     = {[DOI]}
}
```

Code (version DOI): [Zenodo version DOI]
Data: [Zenodo DOI]

## 17. Contact and acknowledgements

Questions and bug reports: please open an issue at https://github.com/LAB/REP/issues
or contact [NAME, EMAIL].

Funding: [FUNDING]. [Other acknowledgements.]
