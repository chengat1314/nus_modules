/* C Example */
#include <stdio.h>
#include <mpi.h>


int main (int argc, char** argv) {
  int rank, size;

  MPI_Init (&argc, &argv); /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */
  double a = 1.0;
  int rankid;
//  for (int i=0; i<1e10; ++i) 
//    a += a*a/3;
  if(rank == 0){
	  
	MPI_Recv(&rankid,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG ,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	printf( "The process 0 receive from process %d\n", rankid);
  }else{
	int rankid = rank;
	MPI_Send(&rankid,1,MPI_INT,0,0,MPI_COMM_WORLD);
	printf( "The process %d of %d sent message to process 0\n", rank, size );
  }
  //printf( "result %f \n", a );
  MPI_Finalize();
  return 0;
}

