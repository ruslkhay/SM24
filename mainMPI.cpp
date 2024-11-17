#include <iostream>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  int numProc;
  int rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::printf("%d: hello (p=%d)\n", rank, numProc);
  std::printf("%d: goodbye\n", rank);

  MPI_Finalize();
  return 0;
}