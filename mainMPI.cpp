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

void Send(Grid::line_t boardVals, int procNum, int rank, int nextRank) {

  std::pair<int, int> buffSize(boardVals.size(), 1);
  MPI_Send(&buffSize, 2, MPI_INT, nextRank, 1, MPI_COMM_WORLD);

  MPI_Send(&boardVals[0], buffSize.first, MPI_DOUBLE, nextRank, 0,
           MPI_COMM_WORLD);
}

void _Receive(int rank, int prevRank) {
  std::pair<int, int> bS;
  MPI_Recv(&bS, 2, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::vector<double> storage(bS.first);
  MPI_Recv(&storage[0], bS.first, MPI_DOUBLE, prevRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  debugReceivePrint(storage, rank, prevRank, bS);
}

Grid::line_t Receive(int rank, int prevRank) {
  std::pair<int, int> bS;
  MPI_Recv(&bS, 2, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  Grid::line_t storage(bS.first);
  MPI_Recv(&storage[0], bS.first, MPI_DOUBLE, prevRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  debugReceivePrint(storage, rank, prevRank, bS);
  return storage;
}

// const int M = 4, N = 4;
const int M = 10, N = 10;
const int maxIter = 1e5;
const double tolerance = 1e-6;
const double h1 = 3.0 / M, h2 = 3.0 / N;
auto method = sMethod::lin;

int main(int argc, char **argv) {
  int rank, size;
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Barrier(MPI_COMM_WORLD);
  std::pair<int, int> buffSize(0, 0);
  int nextRank = (rank + 1) % size;
  int prevRank = rank == 0 ? size - 1 : rank - 1;

  // Communicate processes
  if (rank % 2 == 0) {
    auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
    auto solution =
        Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
    solution.Find(method);
    auto boardVals = solution.GetRightBoarder();
    debugSendPrint(solution.GetNodes(), rank, nextRank, x0, xM, y0, yN,
                   boardVals);
    Send(boardVals, size, rank, nextRank);
    auto boarderVal = Receive(rank, prevRank);
    // if (size == 2) {
    //   solution.SetRightBoarder(boarderVal);
    // }
  } else {
    auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
    printf("(%d; %d), (%d, %d)\n", x0, xM, y0, yN);
    auto solution =
        Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
    auto boarderVal = Receive(rank, prevRank);
    // if (size == 2) {
    //   solution.SetLeftBoarder(boarderVal);
    // }
    solution.Find(method);
    auto boardVals = solution.GetRightBoarder();
    debugSendPrint(solution.GetNodes(), rank, nextRank, x0, xM, y0, yN,
                   boardVals);
    Send(boardVals, size, rank, nextRank);
  }
  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}