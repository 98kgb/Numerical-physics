# -*- coding: utf-8 -*-
"""
Test for the conversion (Matlab to python)

@author: Gibaek
"""

import numpy as np

C = 299792458    # Speed of light [m/s]

save_enable = 1 # Saving options, 1 stands for save, 0 doesn't save. always put it 1
isMat = 0
plot_evo_enable = 1
plot_spectrum_enable = 1
plot_beta2_enable = 1


isTM = 0;
Nsim = 3;  # set simulation iteration times
lambda0 = 1300*1e-9;   #m
Ch = 0;# Initial chirp
betas = [-1.172554609910727e-25, -4.291112339397303e-39];
betas = False;

TaperDispersion; # TaperDispersion script check, for dispersive

P0_list=[600];
