# INME: a canonical neural circuit for motion direction estimation

INME (intrinsically nonlinear motion estimation) is a neural circuit model that
estimates the direction of visual motion using only cells selective to cardinal
directions, with elements that emulate dendritic nonlinearities and
backpropagating action potentials.

This repository contains everything needed to reproduce the results: the model,
the generators for all the stimuli, the scripts that run the eleven experiments,
the code that evaluates five optical-flow networks on the same stimuli, and the
scripts that draw the figures and build the tables.

## Requirements

| | |
|---|---|
| MATLAB | R2021b or newer, with the Image Processing Toolbox |
| Python | 3.9 or newer, with `numpy`, `scipy`, `h5py` and `matplotlib` |
| For the networks | `torch`, `opencv-python`, and the published weights of each network |

Running the model needs only MATLAB. Python is needed to draw the figures and to
evaluate the networks. Nothing else has to be downloaded: every stimulus is
generated from code.

## Layout

| Folder | What it holds |
|---|---|
| `model_stinrf/` | The model (`stINRF3_compute*.m`), its configuration (`stINRF3_config.m`), the readouts, the stimulus generators (`make_*.m`), the experiment scripts (`run_*.m`) and the figure script `plot_figures.py` |
| `stimuli_dnn/` | The same stimuli rendered for the networks. `dnn_v25_config.m` holds the sampling of each experiment, `export_dnn_v25_*.m` write the stimulus sets |
| `dnn_inference/` | The network runners: `run_dnn_fullvideo_gauss.py` for RAFT, PWC-Net and LiteFlowNet2, `run_multiframe_pipeline.py` for ME-Attention and Dual, sharing the readout in `stinrf3_readout.py` |
| `pipeline/` | Scripts that generate a whole stimulus set (`gen_*.m`), pass it through the networks (`run_dnn_experiment.sh`) and collect the model's per-trial responses (`run_*.m`) |
| `analysis/` | `build_tables.py` builds the result tables, `make_videos.py` renders the stimuli to video |

Three helpers sit at the root. `inme_setup.m` puts the repository on the MATLAB
path and returns the output folder; call it before anything else. `inme_output.m`
returns that folder on its own, and `inme_root.m` the repository path.

## Where results are written

Everything goes to `output/` inside the repository. To send it somewhere else,
set `INME_OUT` before you start:

```bash
export INME_OUT=/data/inme
```

The shell scripts call `python3`. If the networks live in a virtual environment,
point to its interpreter with `INME_PYTHON`:

```bash
export INME_PYTHON=~/venv_dnn/bin/python
```

## Tutorial

### 1. Run one experiment

Each experiment is a script in `model_stinrf/`. Run it from the repository root;
`inme_setup` puts the code on the MATLAB path and has to be called first:

```bash
cd INME
matlab -batch "inme_setup; run_ReversePhi_1D"
```

It prints its results, writes them under `output/results/ReversePhi/` and exports
what its figure needs to `output/export/`.

These are the eleven experiments of the paper and the script that runs each one:

| Experiment | Script |
|---|---|
| Missing-fundamental illusion | `run_Missing_1D` |
| Reverse-Phi illusion | `run_ReversePhi_1D` |
| Contrast saturation | `run_Contrast_1D` |
| Motion masking | `run_Masking_nomask_1D`, then `run_Masking_1D` |
| Second-order motion | `run_SecondOrder_1D` |
| Direction of gratings | `run_Gratings` |
| Direction of Type I plaids | `run_PlaidI` |
| Direction of Type II plaids | `run_PlaidII` |
| Interaction among scales | `run_Interaction_1D` |
| Breaking of time-reversal symmetry | `run_LightEdges`, `run_Uneven`, `run_Opponency_1D` (control) |
| Static motion illusions | `run_snakes_2D`, `run_CentralDrift_2D`, `run_RotatingSnakesPeripheral_2D` |

Run the masking control before masking itself: the masked responses are reported
relative to the unmasked one.

### 2. Draw its figure

```bash
python3 model_stinrf/plot_figures.py ReversePhi
```

It finds the export folder on its own, in the same place the experiments wrote
it. Name the figures you want, or name none to redraw every figure that has been
exported. `--export DIR` reads from somewhere else and `--out DIR` writes the
images somewhere else.

### 3. Run everything

```bash
matlab -batch "inme_setup; run_all_experiments"
python3 model_stinrf/plot_figures.py
```

This rebuilds all the model results from nothing. It takes a while; to redo a
single figure, use step 1 instead.

### 4. Run the networks on the same stimuli

First render the stimuli at the sampling each network needs. The networks only
behave well when the displacement between consecutive frames stays between 10 and
40 pixels, so the pixels per degree and the frame rate are chosen per experiment,
while the physical stimulus — degrees per second, cycles per degree, duration —
is the same one the model sees. All of that lives in `dnn_v25_config.m`.

```bash
matlab -batch "inme_setup; gen_reps"
```

Then pass the set through the five networks:

```bash
pipeline/run_dnn_experiment.sh output/stimuli/ReversePhi output/results/ReversePhi
```

Each network writes one direction and one magnitude per stimulus, read out with
the same Gaussian pooling and temporal weighting the model uses, so the two sets
of estimates are directly comparable.

### 5. Build the tables

```bash
python3 analysis/build_tables.py
```

It reads the per-trial results under `output/` and writes the result tables.

## How the model is organised

The model has three stages, and `stINRF3_compute.m` applies them in order:

1. A linear spatio-temporal filter, Gaussian in space and biphasic in time,
   followed by a pointwise nonlinearity.
2. The motion sensor itself, which is where the nonlinearity lives: a Gabor
   filter bank tuned to cardinal directions, combined through a saturating
   function whose argument is the difference between the local signal and a
   pooled version of it.
3. A readout, Gaussian in space and weighted in time, which reduces the response
   to one vector per stimulus.

`stINRF3_config.m` holds the parameters of every experiment in one place, so
changing a stimulus and changing the model are separate operations.
`stINRF3_compute_1D.m` is the one-dimensional version, used by the experiments
whose stimulus varies along a single spatial axis.

## License

Released under the MIT License. See `LICENSE`.
