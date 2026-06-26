# Pulsar Nulling Fraction Calculator

A tool for measuring the **nulling fraction** of radio pulsars from single-pulse data. Nulling is a phenomenon where a pulsar intermittently switches its radio emission off for one or more pulse periods. This code measures what fraction of pulses are "nulls" using two methods:

1. **Wang method** — fits the negative tail of the on-pulse flux distribution to estimate the null fraction
2. **Bayesian (MCMC) method** — models the on-pulse flux distribution as a mixture of a null component (Gaussian noise) and an emission component (log-normal), then samples the posterior using `emcee`

## Requirements

- Python 3.x
- Dependencies listed in `requirements.txt`: `numpy`, `scipy`, `matplotlib`, `emcee`, `corner`

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running with the example data

Example data for pulsar **J1559-5545** is included in `pulsar_data/J1559-5545/`. A ready-made bash script is provided to run the calculator on this data:

```bash
bash bash_script/nf_calculator_J1559-5545.sh
```

This is equivalent to running:

```bash
python python_code/nf_calculator.py \
    -pulsar J1559-5545 \
    -fb 114 \
    -lb 139 \
    -outdir pulsar_results/ \
    -datadir pulsar_data/ \
    -bins 512
```

Results and plots will be written to `pulsar_results/J1559-5545/`.

## Arguments

| Argument   | Description                                           | Example         |
|------------|-------------------------------------------------------|-----------------|
| `-pulsar`  | Pulsar name (must match the subdirectory in datadir)  | `J1559-5545`    |
| `-fb`      | First phase bin of the on-pulse window                | `114`           |
| `-lb`      | Last phase bin of the on-pulse window                 | `139`           |
| `-outdir`  | Directory where output plots and files are saved      | `pulsar_results/` |
| `-datadir` | Directory containing the pulsar data subdirectories   | `pulsar_data/`  |
| `-bins`    | Number of phase bins per pulse period                 | `512`           |

## Input data format

Each observation is a two-column ASCII file (`*.sp.asc`) with a comment header on the first line. The columns are phase bin number and flux density:

```
# 51844.0 15961.0000001 0.9573023930 1 1516.500 212.900 512 7 1 J1559-5545
1 -0.111111
2 -2.428571
...
```

Data files should be placed in `pulsar_data/<pulsar_name>/`.

## Reproducibility

The Wang method is fully deterministic. The Bayesian method uses MCMC (`emcee`) and is stochastic: walker starting positions are randomly perturbed and the sampler draws different chains each run, so the Bayesian nulling fraction and its uncertainties will vary slightly between runs. With the default settings (100 burn-in steps, 1000 sample steps, 10 walkers) this variation is typically small, but to get exactly reproducible results you can add `np.random.seed(<integer>)` near the top of `nf_calculator.py` before the MCMC section.

## Outputs

For each observation file, the following plots are saved:

| File | Description |
|------|-------------|
| `*_waterfall_mean.png` | Waterfall plot of single pulses and the mean profile |
| `*_flux.png` | Per-pulse on-window flux density time series |
| `*_hist_before.png` | On- and off-pulse flux histograms before scaling |
| `*_hist_neg.png` | Histogram of the negative flux values used by the Wang method |
| `*_hist_scaled.png` | Scaled histogram showing Wang method fit |
| `*_bayes_fit.png` | Histogram with Bayesian model overlay |
| `*_corner_plot.png` | MCMC posterior corner plot |

Summary results across all observations are saved to:

**`<outdir>/<pulsar>/<pulsar>.txt`** — one row per quantity, one column per observation (sorted chronologically by MJD):

| Row | Quantity | Description |
|-----|----------|-------------|
| 1 | MJD | Modified Julian Date of the observation |
| 2 | Bayesian NF | Bayesian nulling fraction (posterior median) |
| 3 | NF lower uncertainty | Lower uncertainty on the Bayesian NF |
| 4 | NF upper uncertainty | Upper uncertainty on the Bayesian NF |
| 5 | S/N proxy | Signal-to-noise proxy for the observation |
| 6 | Number of profiles | Number of single pulses in the observation |
| 7 | Wang NF | Wang method nulling fraction |
| 8 | Max consecutive nulls | Longest consecutive null run in the observation |

**`<outdir>/<pulsar>/<pulsar>_consecutive_nulls_and_non.txt`** — null and non-null run lengths
