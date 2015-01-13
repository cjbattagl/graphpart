make
mpirun -np 9 generator_test_mpi 7
cat file* > mygraph.mat
rm file*
