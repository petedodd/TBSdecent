#!/bin/bash

# CEA with restricted costs
R --slave --vanilla --args <decent_run.R CEA DECENT_EQUIPMENT
R --slave --vanilla --args <decent_run.R CEA DECENT_PERSONNEL
R --slave --vanilla --args <decent_run.R CEA DECENT_SUPPLY
R --slave --vanilla --args <decent_run.R CEA DECENT_TRAINING


# main CEA
# NOTE probably makes sense to run after variants in case any variant outputs are not labelled as such
R --slave --vanilla --args <decent_run.R CEA


# main BIA
R --slave --vanilla --args <decent_run.R BIA

# NOTE
# if adding for BIA:
# R --slave --vanilla --args <decent_run.R BIA DECENT_TRAINING
# with inputs files like TBS.DECENT_TRAINING.costsBIA.csv etc
