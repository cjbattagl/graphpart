make
#rm -f mygraph*
#rm out_pcsr*
#rm out_tup*
#rm out_csr*
#rm out_permed*
mpirun -np 8 ./graph500_part_test 10 16
mpirun -np 8 ./graph500_part_test 11 16
mpirun -np 8 ./graph500_part_test 12 16
mpirun -np 8 ./graph500_part_test 13 16
mpirun -np 8 ./graph500_part_test 14 16
mpirun -np 8 ./graph500_part_test 15 16
mpirun -np 8 ./graph500_part_test 16 16
mpirun -np 8 ./graph500_part_test 17 16
mpirun -np 8 ./graph500_part_test 18 16
mpirun -np 8 ./graph500_part_test 19 16

#cat out_tup* > mygraphtup.mat
#cat out_csr* > mygraphcsr.mat
#cat out_permed* > mygraphperm.mat
#cat out_pcsr* > mygraphpcsr.mat

