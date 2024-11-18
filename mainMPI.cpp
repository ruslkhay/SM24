#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>

int main(int argc, char *argv[]) {
  int numProc;
  int rank;
  std::vector<int> ar{};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::printf("%d: hello (p=%d)\n", rank, numProc);
  if (rank % 2 == 0) {
    MPI_Send(&rank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
  }
  std::printf("%d: goodbye\n", rank);

  MPI_Finalize();
  return 0;
}