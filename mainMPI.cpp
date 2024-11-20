#include "src/mpi.hpp"
#include "src/solution.hpp"
#include <algorithm>
#include <cstdlib>
#include <mpi.h>
#include <stdio.h>
#include <vector>

/*
1. Give each process corresponding sub-areas. Therefore implement algorithm
2. Calculate solution in sub-area
3. Send corresponding to the process id boarder values
4. Calculate
*/

void Send(std::vector<std::vector<int>> &grid, int procNum, int rank,
          int nextRank, int M, int N) {
  auto [x0, xM, y0, yN] = GetSectors(procNum, rank, M, N);
  debugSendPrint(grid, rank, nextRank, x0, xM, y0, yN);
  auto tmpBuf = Grid(prepareSubGrid(grid, x0, xM, y0, yN));
  std::pair<int, int> buffSize(tmpBuf.GetM(), 1);
  MPI_Send(&buffSize, 2, MPI_INT, nextRank, 1, MPI_COMM_WORLD);
  MPI_Send(&(tmpBuf.GetRightBoarder())[0], buffSize.first, MPI_DOUBLE, nextRank,
           0, MPI_COMM_WORLD);
}

void Receive(int rank, int prevRank, int M, int N) {
  std::pair<int, int> bS;
  MPI_Recv(&bS, 2, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::vector<double> storage;
  storage.reserve(bS.first);
  MPI_Recv(&storage[0], bS.first, MPI_DOUBLE, prevRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  debugReceivePrint(storage, rank, prevRank, bS);
}

const int M = 5, N = 5;
std::vector<std::vector<int>> grid = {{1, 1, -1, 1, 9},
                                      {2, 2, -2, 2, 8},
                                      {3, 3, -3, 3, 7},
                                      {4, 4, -4, 4, 6},
                                      {5, 5, -5, 5, 0}};

int main(int argc, char **argv) {
  int rank, size;
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Fill data for grid only in process 0
  if (rank == size - 1) {
    // Display original grid values
    printf("\nOriginal grid values:\n");
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        printf("%4d", grid[i][j]);
      }
      printf("\n");
    }
    std::cout << "\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);
  std::pair<int, int> buffSize(0, 0);
  int nextRank = (rank + 1) % size;
  int prevRank = rank == 0 ? size - 1 : rank - 1;
  // Communicate processes
  if (rank % 2 == 0) {
    Send(grid, size, rank, nextRank, M, N);
    Receive(rank, prevRank, M, N);
  } else {
    Receive(rank, prevRank, M, N);
    Send(grid, size, rank, nextRank, M, N);
  }
  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}