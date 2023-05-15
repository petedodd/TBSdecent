#!/bin/bash

# CEA with restricted costs
R --slave --vanilla --args <decent_run.R CEA DECENT_EQUIPMENT &
R --slave --vanilla --args <decent_run.R CEA DECENT_PERSONNEL &
R --slave --vanilla --args <decent_run.R CEA DECENT_SUPPLY &
R --slave --vanilla --args <decent_run.R CEA DECENT_TRAINING &


# main CEA
# NOTE probably makes sense to run after variants in case any variant outputs are not labelled as such
R --slave --vanilla --args <decent_run.R CEA

# example SA fixing prevalence
# R --slave --vanilla --args <decent_run.R CEA DECENT 500e-5
# version fixing prevalence and using 6% discount rate
# R --slave --vanilla --args <decent_run.R CEA DECENT 500e-5 0.06
# (output to above prepended: 0.0050.06)


# # set to run the sensitivity analysis grid
# R --slave --vanilla --args <decent_run.R CEA DECENT 500e-5 0.06
# R --slave --vanilla --args <decent_run.R CEA DECENT 500e-5 0.03
# R --slave --vanilla --args <decent_run.R CEA DECENT 500e-5 0.00
# R --slave --vanilla --args <decent_run.R CEA DECENT 200e-5 0.06
# R --slave --vanilla --args <decent_run.R CEA DECENT 200e-5 0.03
# R --slave --vanilla --args <decent_run.R CEA DECENT 200e-5 0.00
# R --slave --vanilla --args <decent_run.R CEA DECENT 50e-5 0.06
# R --slave --vanilla --args <decent_run.R CEA DECENT 50e-5 0.03
# R --slave --vanilla --args <decent_run.R CEA DECENT 50e-5 0.00



# main BIA
# R --slave --vanilla --args <decent_run.R BIA

# NOTE
# if adding for BIA:
# R --slave --vanilla --args <decent_run.R BIA DECENT_TRAINING
# with inputs files like TBS.DECENT_TRAINING.costsBIA.csv etc
