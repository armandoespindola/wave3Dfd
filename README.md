# Wave3Dfd-4TH
Wave3Dfd simulates forward wave propagation 4th-order using MPI-OpenMP

How to run the example:

1) make clean

2) make all

3) copy bin/*.out to Example/

4 cd Example/

4) mpirun -np 24 ./wave3Dfd.out -nFile parameter.par
