#!/bin/bash

# to plot the data vs simulations for IP8 colliding / not colliding bunches
# two data files needs to be in RESULTS_COMPARE_NOMINAL
python plotVerHorOffsetsComparisonWithData.py RESULTS/train_4440_170m_new_fromMatlab.in/ RESULTS_COMPARE_NOMINAL/train_4440_170m.in.only15/ &
