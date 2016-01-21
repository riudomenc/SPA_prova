#!/bin/bash

# 1) run SPA
./spa Prades_ndf.config

# 2) plot output
R --vanilla --slave --args "figures/" FALSE FALSE "ndf" <./plot.R

# argument 5: output figure directory
# argument 6: SDE?
# argument 7: plot parameter correlations (time consuming)?
# argument 8: defoliated or nondefoliated?

