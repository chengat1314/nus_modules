/* C Example */
#include <stdio.h>
#include <mpi.h>


int main (int argc, char** argv) {
  int rank, size;

  MPI_Init (&argc, &argv); /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */
  double a = 1.0;
//  for (int i=0; i<1e10; ++i) 
//    a += a*a/3;
     
  printf( "Hello world from process %d of %d\n", rank, size );
  //printf( "result %f \n", a );
  MPI_Finalize();
  return 0;
}

