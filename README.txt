This repository contains all code to reproduce the tables and figures in "Inference conditional on selection: a review", by Neufeld, Perry, and Witten (2026). 

The code for Figure 2 will take ~20 minutes to run, but can be run on a local machine. 
The code for Figure 1 is intended to be submitted to run on a cluster using submit_array.sh. 

In the Figure 1 code, the files in the Utils folder and CART.py are taken directly from https://github.com/yiling-h/SI-CART/tree/881287397c013a8997f93567606c6457d7a6f547/Simulations/Replication, as they are meant to reproduce the results of Bakshi, Huang, Panigrahi, and Dempsey (2025). The file RRT_fixed.py is also directly inspired by code from that repository. Some code for sample splitting, Full CSI, and the classical method is also re-used from Neufeld, Gao, and Witten (2022), and was taken from https://github.com/anna-neufeld/treevalues-simulations. 