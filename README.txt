This repository contains all code to reproduce the tables and figures in "Inference conditional on selection: a review", by Neufeld, Perry, and Witten (2026). 

The code for Tables 1 and 2 can be run locally, and will run quickly. 

To reproduce Figure 1 from existing results, you can use the file "analyze_results.R", which will draw on results stored in the folder "mar30". To reproduce the results in the folder "mar30", you will need to run the file "run_chunk.R" 100 times. This is quite slow, and it intended to be run on a cluster using submit_array.sh. 

In the Figure 1 folder, the files in the Utils folder and CART.py are taken directly from https://github.com/yiling-h/SI-CART/tree/881287397c013a8997f93567606c6457d7a6f547/Simulations/Replication, as they are meant to reproduce the results of Bakshi, Huang, Panigrahi, and Dempsey (2025). The file RRT_fixed.py is also directly inspired by code from that repository. Some code for sample splitting, Full CSI, and the classical method in one_rep.R is also re-used from Neufeld, Gao, and Witten (2022), and was taken from https://github.com/anna-neufeld/treevalues-simulations. 

To reproduce Figure 2, the code can be run locally in 20-30 minutes. 
