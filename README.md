# Intrinsic Disorder TCR Model

Ideal chain simulation of T Cell Receptor chains, from Clemens, Dushek and Allard "Intrinsic disorder in the T cell receptor creates cooperativity and controls ZAP70 binding" [bioRxiv: 10.1101/2020.05.21.108662v1](https://www.biorxiv.org/content/10.1101/2020.05.21.108662v1).

## Quickstart

Requires C compiler `gcc` and `make` to run the simulation, and Matlab to analyze results.

* In `polymer_sim/`

```
make
./batchMetropolis.sh
```

* This script runs 64 simulations in parallel and takes about 2 hours on a 1.6 GHz Dual-Core Intel Core i5.

* In `analysis/`, run Matlab m-file `Plot_UnweightedAverageBindingRateVSNumBind.m`

This generates Figures S1 and S2 from Clemens et al.,
