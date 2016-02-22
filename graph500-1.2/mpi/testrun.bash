make
mpirun -np 8 ./graph500_part_test 12 16
#mpirun -np 32 ./graph500_mpi_simple 16 16
cat out_tup* > mygraphtup.mat
cat out_csr* > mygraphcsr.mat
cat out_permed* > mygraphperm.mat
cat out_pcsr* > mygraphpcsr.mat
rm out_pcsr*
rm out_tup*
rm out_csr*
rm out_permed*
