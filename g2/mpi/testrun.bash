make
#rm -f mygraph*
#rm out_pcsr*
#rm out_tup*
#rm out_csr*
#rm out_permed*
mpirun -np 16 ./graph500_mpi_simple 15 16
mpirun -np 16 ./graph500_mpi_simple 16 16
mpirun -np 16 ./graph500_mpi_simple 17 16
mpirun -np 16 ./graph500_mpi_simple 18 16
mpirun -np 16 ./graph500_mpi_simple 19 16

#cat out_tup* > mygraphtup.mat
#cat out_csr* > mygraphcsr.mat
#cat out_permed* > mygraphperm.mat
#cat out_pcsr* > mygraphpcsr.mat

