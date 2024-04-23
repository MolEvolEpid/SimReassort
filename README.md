# SimReassort
Code from Lin et al., "The number and pattern of viral genomic reassortments are not necessarily identifiable from segment trees", Mol Biol Evol, 2024; msae078. DOI: [https://doi.org/10.1093/molbev/msae078](https://doi.org/10.1093/molbev/msae078).

**SimReassort** is a minimal package modified from [**phylopomp**](https://github.com/kingaa/phylopomp), which simulates and infers from genealogies given a broad classes of compartmental epidemiological models, described in King et al (2022), doi: [https://doi.org/10.1016/j.tpb.2021.11.003](https://doi.org/10.1016/j.tpb.2021.11.003).

**SimReassort** simulates transmission dynamics and the genealogical evolutions for two-segment viruses, using a linear-birth-death model structure.

## Directory structure

-   `R`: contains `R` codes for simulating and ploting genealogies, dismissing invisible reassortments, and computing the remove-and-rejoin table

-   `examples`: contains examples for simulating genealogies and computing remove-and-rejoin table

-   `src`: contains `C/C++` codes modified from package `phylopomp`

-   `man`: contains manuals for available functions

-   `data`: contains simulated data

## Installation

A list of packages are required: `ggtree`, `ggplot2`, `ape`, `dplyr`, `tibble`, `foreach`, `doParallel`, `tidyr`, `pomp`, `cowplot`, `scales`, and `stringr`.

``` r
install.packages("devtools")
devtools::install_github("MolEvolEpid/SimReassort")
```
