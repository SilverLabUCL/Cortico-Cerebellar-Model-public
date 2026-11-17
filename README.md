# [PAPER TITLE]: simulation and analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

C++ simulation code and MATLAB analysis scripts accompanying the manuscript
**"[PAPER TITLE]"** by [AUTHORS] ([preprint link / DOI], submitted to [JOURNAL]).

> **TODO before publishing:** replace every `[PLACEHOLDER]` and `TODO` in this file.
> Items marked **[CONFIRM]** are my best reading of your notes and should be checked.

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
- text files listing the exact commands used for every simulation and figure in the paper,
- a Dockerfile for a reproducible build environment.

The repository does **not** contain the raw simulation output. The processed data
behind each figure are archived on Zenodo (see [Data](#10-data-zenodo)).

## 2. Repository structure

```
REP/
├── RNN_only.cpp            # simulation: RNN only
├── RNN+FFL.cpp             # simulation: RNN + FFL
├── headers/                # C++ header files
├── simulation_results/     # main output files (created/filled by the programs)
├── connectivity/           # network weights; only written for single example simulations
├── *.m                     # MATLAB scripts (loading, statistics, plotting)
├── [commands_*.txt]        # TODO: commands for all paper simulations and figures
├── Dockerfile
├── LICENSE
└── README.md
```

TODO: if you decide to move the MATLAB scripts into their own folder (for example `matlab/`),
update this tree and [Section 8](#8-matlab-workflow). Keeping them next to the C++ code, as now, is fine.

## 3. Requirements

### Tested environment

| Component | Version |
|---|---|
| Operating system | Ubuntu 24.04 (Windows is not supported) |
| Compiler | g++ 13.3.0 |
| [Armadillo](https://arma.sourceforge.net/) | 14.4.0 |
| [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) | 2.8 |
| [Boost](https://www.boost.org/) | 1.83.0 (components used: `program_options`, `version`) |
| BLAS / LAPACK | pthread-enabled OpenBLAS (used for both, accessed through Armadillo) |
| OpenMP | provided by g++ (`-fopenmp`), used by `RNN+FFL` **[CONFIRM: also RNN_only?]** |
| [nlohmann/json](https://github.com/nlohmann/json) | [VERSION] **[CONFIRM: bundled in `headers/`?]** |
| MATLAB | R2024b (24.2.0.3212159), Update 9 |
| MATLAB toolboxes | Statistics and Machine Learning Toolbox |

MATLAB scripts read the JSON output with `jsondecode`.

### Installing the dependencies

On Ubuntu 24.04 the following packages provide the compiler, Boost and OpenBLAS:

```bash
sudo apt update
sudo apt install build-essential libboost-program-options-dev libopenblas-dev liblapack-dev
```

> **TODO:** to my knowledge Ubuntu 24.04's packaged versions of Armadillo and GSL are older than
> the ones listed above (Armadillo 12.x, GSL 2.7.x). If you built these from source, add the
> instructions here, otherwise `apt install libarmadillo-dev libgsl-dev` will give slightly
> different versions. The [Docker image](#13-docker) pins the exact versions.

## 4. Compilation

Compile from the repository root:

```bash
g++ -O3 -o RNN_only.out RNN_only.cpp -lgsl -lgslcblas -lm -larmadillo -lboost_program_options

g++ -O3 -o RNN+FFL.out RNN+FFL.cpp -lgsl -lgslcblas -lm -larmadillo -lboost_program_options -fopenmp
```

| Flag | Purpose |
|---|---|
| `-O3` | highest standard optimisation level (important for simulation speed) |
| `-lgsl -lgslcblas -lm` | link GSL, its CBLAS interface and the math library |
| `-larmadillo` | link Armadillo (which uses OpenBLAS/LAPACK underneath) |
| `-lboost_program_options` | link Boost.Program_options (command-line parsing) |
| `-fopenmp` | enable OpenMP multithreading |

Expected compilation time: [TODO: seconds/minutes per program].

Before running `RNN+FFL.out`, set the number of OpenMP threads (the paper simulations used 6):

```bash
export OMP_NUM_THREADS=6
```

## 5. Running the program

Both programs take their settings from the command line, parsed with
[Boost.Program_options](https://www.boost.org/doc/libs/release/doc/html/program_options.html).
The options used in the test simulations are:

| Option | Meaning | Values used here |
|---|---|---|
| `--learnmode` | what the program does | `none` (single simulation), `scan_lambda` (parameter sweep over `lambda`) **[CONFIRM: list all modes]** |
| `--Ncomp` | [TODO: meaning] | `3` (single simulation), `1` (sweep) |
| `--Nreal` | number of network realisations (repetitions) | `1`, `3`; the paper uses 50 |
| `--seed` | seed of the random number generator | `-1` **[CONFIRM: what does -1 mean?]** |
| `--lambda` | [TODO: meaning] | `10` |
| `--Floquet` | compute Floquet multipliers (`ON`/`OFF`) | `OFF` |

The full list of options, their defaults and descriptions is printed by:

```bash
./RNN_only.out --help
./RNN+FFL.out --help
```

> **TODO (parameter descriptions):** I'd recommend keeping `--help` as the single source of truth,
> keeping only the essential options in this table, and putting longer descriptions
> (including the contents of `parameter.dat` and `flag.dat`) in a separate
> `docs/PARAMETERS.md` linked from here.

Every run also writes a **log file** ([TODO: file name and location]) that accumulates a
record of which stage of the program is currently executing (for example `Fourier analysis`).
This is useful for tracking progress in long simulations.

## 6. Test simulations

Run these from the repository root after compiling. All four commands print a short
console summary: a few key parameters, the names of the folders where output is saved,
and the Armadillo and Boost versions in use.

| # | Command | Runtime | What it produces |
|---|---|---|---|
| 1 | `./RNN_only.out --learnmode none --Ncomp 3 --Nreal 1 --seed -1 --lambda 10 --Floquet OFF` | ~5 s | A single RNN-only simulation like the one in **Fig. 3b-d** |
| 2 | `./RNN+FFL.out --learnmode none --Ncomp 3 --Nreal 1 --seed -1 --lambda 10 --Floquet OFF` | ~5 s | A single RNN+FFL simulation like the one in **Fig. 3f-h** |
| 3 | `./RNN_only.out --learnmode scan_lambda --Ncomp 1 --Nreal 3 --seed -1 --Floquet OFF` | ~6 min | RNN-only parameter sweep: a low-repetition version of **Fig. 1d** |
| 4 | `./RNN+FFL.out --learnmode scan_lambda --Ncomp 1 --Nreal 3 --seed -1 --Floquet OFF` | ~8 min | RNN+FFL parameter sweep: a low-repetition version of **Fig. 1d** |

Runtimes are for the hardware described in [Section 11](#11-computational-resources-and-runtime).

### What the output should look like

- **Tests 1 and 2** (single simulations) write to both `simulation_results/` and
  `connectivity/`. See [Output files](#7-output-files) for the contents.
- **Tests 3 and 4** (sweeps) write only to `simulation_results/`. For each parameter set `i`
  of the sweep, a folder `parameter_i/` is created containing the same files as in the
  single-simulation case.

### Checking that a test worked

1. The console summary appears and reports the expected Armadillo (14.4.0) and Boost (1.83.0) versions.
2. The folders and files listed in [Output files](#7-output-files) exist and are non-empty.
3. Load the result in MATLAB ([Section 8](#8-matlab-workflow)) and run [TODO: name of plotting script]
   to compare with Fig. 3b-d (tests 1) or Fig. 3f-h (test 2).

> **TODO:** consider adding a small reference output and a check script (with a tolerance) so that
> users can verify their installation automatically. Also note that with `--seed -1` the
> run may not be reproducible (see [Section 12](#12-reproducibility-notes)).

## 7. Output files

Most simulation results are stored as **JSON** (written with
[nlohmann/json](https://github.com/nlohmann/json)); [TODO: complete: which files are JSON and which are plain text].

**`simulation_results/`**

| File / folder | Contents |
|---|---|
| `parameter.dat` | all simulation parameters |
| `flag.dat` | all flags |
| [TODO: file names] | population firing rates from the test phase and from an open-loop phase |
| `parameter_1/` (`parameter_i/` in sweeps) | for each network realisation (`Nreal`): test-phase readout neuron activity, target signals, network firing-rate statistics, network representational measures, and example single-neuron time constants from a single trial |

**`connectivity/`** (only for single example simulations): network weight files (recurrent
weights, readout weights, etc.).

Each new run [TODO: overwrites / appends to / creates a new subfolder in] these folders.

## 8. MATLAB workflow

```
C++ programs  ->  JSON/text files  ->  MATLAB loading (jsondecode -> structs)
              ->  statistics over simulations  ->  figures
```

The MATLAB scripts are located in the same folder as the C++ programs. Run them with the
repository root as the current folder (or add it with `addpath`).
Requires MATLAB R2024b and the Statistics and Machine Learning Toolbox.

| Script | Purpose | Input | Output |
|---|---|---|---|
| [`load_results.m`] | reads the JSON files with `jsondecode` and converts them into MATLAB structs | `simulation_results/` | structs in the workspace |
| [`compute_statistics.m`] | statistics over simulations/realisations | structs | [TODO] |
| [`plot_figX.m`] | generates the figure panels | structs / Zenodo data | figure panels |

TODO: replace the placeholder rows with the real script names (one line per script).

## 9. Reproducing the paper

The exact commands used for all simulations and figures in the paper are listed in text files
included in this repository:

| Figure | Command file | Data (Zenodo) | MATLAB script |
|---|---|---|---|
| Fig. 1 | [`commands_fig1.txt`] | [file] | [script] |
| Fig. 3 | [`commands_fig3.txt`] | [file] | [script] |
| ... | ... | ... | ... |
| Fig. 8 | [`commands_fig8.txt`] | [file] | [script] |

There are two ways to proceed:

- **(a) Regenerate the figures from the processed data (fast).** Download the data from Zenodo
  ([Section 10](#10-data-zenodo)) and run the plotting scripts. No C++ compilation is needed.
- **(b) Rerun the simulations (slow).** Compile the programs, run the commands in the relevant
  command file, then run the MATLAB scripts on the new output. Expect runtimes of hours
  per figure ([Section 11](#11-computational-resources-and-runtime)).

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

All simulations were run on a Lenovo ThinkPad with an Intel Core Ultra 7 processor
([TODO: exact model, RAM, OS]), using 6 OpenMP threads. No cluster was used.

| Task | Approximate runtime |
|---|---|
| Test simulations 1 and 2 | ~5 s each |
| Test simulation 3 | ~6 min |
| Test simulation 4 | ~8 min |
| Full parameter sweep for most figures (50 repetitions, Floquet multipliers computed), per network architecture | [TODO: "a few hours"?] |
| Fig. 8 (most complex architecture, Fig. 8g) | ~50 hours |

Memory: [TODO]. Disk space for the raw output of a full run: [TODO].

## 12. Reproducibility notes

- Results are deterministic given a random number seed (`--seed`).
- [TODO: what does `--seed -1` do? If it draws a random seed, tell users to pass a fixed positive
  seed when they need identical output, and say how the seed is recorded in `parameter.dat`.]
- [TODO: does the number of OpenMP threads affect the numerical results? Small differences
  across machines/BLAS versions are possible if so.]

## 13. Docker

A Dockerfile with all dependencies pinned to the versions above is provided.

```bash
docker build -t rep .
docker run --rm -it -v "$(pwd)":/work rep
```

TODO: adapt these commands to the final Dockerfile, and state whether MATLAB is included
(it is not, since it is proprietary software).

## 14. Troubleshooting

- **`cannot find -lgslcblas` / `-lgsl`:** GSL is missing or not on the linker path.
- **`cannot find -lboost_program_options`:** install `libboost-program-options-dev`.
- **Armadillo link or LAPACK errors:** make sure OpenBLAS/LAPACK development packages are installed.
- **Single-threaded runs of `RNN+FFL.out`:** the program was probably compiled without `-fopenmp`, or `OMP_NUM_THREADS` is not set.
- **Files or folders not found:** run the programs from the repository root so that `simulation_results/` and `connectivity/` resolve correctly [TODO: confirm].

## 15. License

- **Code:** MIT (see [`LICENSE`](LICENSE)).
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

Code (version DOI): [TODO: Zenodo version DOI]
Data: [TODO: Zenodo DOI]

TODO: add a `CITATION.cff` file so that GitHub displays a "Cite this repository" button.

## 17. Contact and acknowledgements

Questions and bug reports: please open an issue at https://github.com/LAB/REP/issues
or contact [NAME, EMAIL].

Funding: [TODO]. [Other acknowledgements.]
