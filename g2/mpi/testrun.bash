make
rm -f mygraph*
rm out_pcsr*
rm out_tup*
rm out_csr*
rm out_permed*

mpirun -np 8 ./graph500_part_test 22 16
mpirun -np 8 ./graph500_mpi_simple 22 16
cat out_tup* > mygraphtup.mat
cat out_csr* > mygraphcsr.mat
cat out_permed* > mygraphperm.mat
cat out_pcsr* > mygraphpcsr.mat

