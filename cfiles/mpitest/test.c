#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void do_stupid_gather();
void do_gather();


int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  //do_stupid_gather();
  do_gather();
}


// do an all gather using the MPI command
void do_gather() {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);
  int localsize = 4;
  
  int* data = (int*)malloc(sizeof(int) * numprocs * localsize);
  int* localdata = (int*)malloc(sizeof(int) * localsize);
  
  for (int i=0; i<numprocs*localsize; i++) { data[i] = 0; }
  for (int i=0; i<localsize; i++) { localdata[i] = rank; }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Send half of localdata. Collective should place it in data in a strided manner
  // (important for dist partitioner)
  int offset = localsize/2;
  
  MPI_Allgather(localdata+offset, offset, MPI_INT, data+offset, localsize, MPI_INT, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Process %d on %s out of %d, ", rank, processor_name, numprocs);
  for (int i=0; i<numprocs * localsize; i++) { printf("%d ", data[i]); }
  printf("\n");
}

// do a stupid allgather by broadcasting from every process using an offset
void do_stupid_gather() {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);
  int* data = (int*)malloc(sizeof(int) * numprocs);
  for (int i=0; i<numprocs; i++) { data[i] = 0; }
  data[rank] = rank;
  printf("Process %d on %s out of %d, ", rank, processor_name, numprocs);
  for (int i=0; i<numprocs; i++) { printf("%d ", data[i]); }
  printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i=0; i<numprocs; i++) { MPI_Bcast(data+i,1,MPI_INT,i,MPI_COMM_WORLD); }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Process %d on %s out of %d, ", rank, processor_name, numprocs);
  for (int i=0; i<numprocs; i++) { printf("%d ", data[i]); }
  printf("\n");
}