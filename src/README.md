A C++ library which supposedly optimizes the force field constant for CHARMM model, regarding the new molecules.

Written by Nguyen Doan Dai, under supervision of Dr. Alexey Alexandrov. (c) 2021.

# Build and test

Assuming one is in the directory ./Summer internship 2021

1. run `cmake --build build`. There should be a file named CHARMM_Optimizer.exe in Windows, for instance.

2. run `cd build` and `./CHARMM_Optimizer` on Linux (or `.exe` on Windows), followed by the url pointing to data file. Ex: one has the file
   ./Summer internship 2021/ethane_dihe_c1_c2.dat, then one runs
   `./CHARMM_Optimizer ethane_dihe_c1_c2.dat`

    See `-h` for more information.

3. The result will be printed into the file log.txt, containing the final RMSE error, force constants of each angles, corresponding to harmonics as defined in multiplicities[], the average, the standard deviation, and the 95%-confidence of the run. Also, for analysis purpose, runtime will also be printed.
