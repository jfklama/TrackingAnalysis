# TrackingAnalysis

Analysis of ILD tracking performance for Long-Lived Particle (LLP) detection at the ILC.

## Overview

Long-Lived Particles travel a macroscopic distance before decaying, leaving a displaced secondary vertex inside the detector. Reconstructing these vertices efficiently and with low fake rate is a key challenge for BSM searches. This repository contains the full analysis chain — from event simulation to efficiency histograms and exclusion limits — for several LLP models studied in the context of the ILD detector at the ILC.

## Physics models

| Model | Abbreviation | Kinematic regime | Decay signature |
|---|---|---|---|
| Inert Doublet Model | IDM | ΔM = 1-5 GeV | A → Z* l+l− |
| Feebly-Interacting Massive Particle | FIMP | M = 60–110 GeV | χ → l+l−ν |
| Axion-Like Particle | ALP | M = 0.3–10 GeV | a → γγ |
| Two-Real-Singlet Model | TRSM | M = 0.4–60 GeV | S → bb̄, ττ |

## Software stack

- **DD4hep / ddsim** — Geant4-based detector simulation
- **Marlin** — modular event reconstruction framework
- **pyLCIO** — Python bindings for the LCIO event data model
- **ROOT / PyROOT** — histogramming and plotting
- **numpy / scipy** — numerical utilities
- **pytest** — unit tests

The full ILD software environment is necessary to run this code. It is available via [iLCSoft](https://github.com/iLCSoft).

## Analysis pipeline

```
1. Generation     WHIZARD / particle gun  →  .stdhep / HEPEvt files
2. Simulation     ddsim (ddsim_steer.py)  →  .slcio with MC truth
3. Reconstruction Marlin (MarlinStdReco)  →  tracks, vertices, PFOs
4. LLP refit      LLPRefitProcessor       →  secondary vertex candidates
5. Analysis       run_*_analysis.py       →  ROOT histograms, efficiency tables
6. Limits         plot_limits.py       →  95% CL exclusion contours
```

## Repository structure

| Path | Role |
|---|---|
| `tracking.py` | Track-to-MC matching, efficiency histograms |
| `vertex_finding.py` | Secondary vertex matching, V0/conversion veto |
| `idm_analysis.py` | Main analysis loop for the IDM model |
| `overlay_analysis.py` | Main analysis loop for the SM backgrounds |
| `sm_analysis.py`, `v0_analysis.py` | Legacy Standard Model background and V0 validation |
| `run_idm_analysis.py` etc. | Batch launchers for each physics model |
| `utils.py` | Geometry helpers, helix calculations, invariant mass |
| `constants.py` | Detector geometry and particle mass constants |
| `histo_class.py` | ROOT histogram container |
| `counters.py` | Event statistics and weighted-event accumulators |
| `ddsim_steer.py` | DD4hep simulation steering file |
| `ParticleLists/` | Custom particle tables for each model/mass point |
| `LLPFinder/` | C++ Marlin processor — LLP vertex finder |
| `LLPRefitProcessor/` | C++ Marlin processor — track refit to secondary vertex |
| `simulate_idm.sh` etc. | Simulation batch scripts |
| `generate_kaons.sh`, `generate_lambdas.sh` | V0 sample generation |
| `setup.sh` | Environment setup |
| `tests/` | Python unit tests |

## Usage

### Simulation

Particle gun (e.g. lambdas):

```shell
ddsim --steeringFile ddsim_steer.py \
      --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
      --enableGun --gun.particle lambda --gun.energy 20*GeV \
      --gun.distribution uniform --gun.isotrop True \
      -N 1000 --outputFile output.slcio
```

From pre-generated input samples:

```shell
ddsim \
    --inputFiles data/idm_test_generation_dM_10.stdhep \
    -N 10000 \
    --compactFile $lcgeo_DIR/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml \
    --physics.pdgfile ParticleLists/particle_idm_dM_10.tbl \
    --outputFile data/idm_test_simulation_dM_10.slcio
```

### Reconstruction

```shell
Marlin MarlinStdReco.xml \
    --constant.lcgeo_DIR=$lcgeo_DIR \
    --constant.DetectorModel=ILD_l5_o1_v02 \
    --constant.OutputBaseName=<outputBaseName> \
    --global.LCIOInputFiles=<inputFileName>.slcio
```

### Analysis

```shell
python run_idm_analysis.py
```

Individual mass points can also be run directly via `idm_analysis.py` by passing the input `.slcio` file path.

## Tests

```shell
pip install pytest
pytest
```

The test suite covers the pure-Python utility functions (geometry, helix math, angular distances) and the event counter class. pyLCIO is mocked so tests run without the full ILD software stack.

## C++ processors

**LLPFinder** identifies displaced vertex topologies by looking for pairs of tracks with large impact parameters and small distances of closest approach. **LLPRefitProcessor** can be used to refit the tracks, potentially improving the efficiency of the **LLPFinder**. Both processors are built with CMake against a Marlin installation. See [LLPRefitProcessor/README.md](LLPRefitProcessor/README.md) for build instructions.
