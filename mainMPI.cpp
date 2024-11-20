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

template <std::size_t M, std::size_t N>
int *prepareSubGrid(std::array<std::array<int, M>, N> &grid, int x0, int xM,
                    int y0, int yN) {
  int count = (xM - x0) * (yN - y0);
  int *tmpBuf = new int[count];
  for (int i = y0; i < yN; ++i) {
    for (int j = x0; j < xM; ++j) {
      tmpBuf[(i - y0) * (xM - x0) + (j - x0)] = grid[i][j];
    }
  }
  return tmpBuf;
}

std::array<int, 4> GetLimits(int rank, int M, int N) {
  const int xMiddle = M / 2 + 1, yMiddle = N / 2 + 1;
  int x0 = 0, xM = 0, y0 = 0, yN = 0;
  switch (rank) {
  case 0:
    x0 = 0, xM = xMiddle, y0 = 0, yN = N - yMiddle;
    break;
  case 1:
    x0 = 0, xM = xMiddle, y0 = N - yMiddle, yN = N;
    break;
  case 2:
    x0 = xMiddle, xM = M, y0 = N - yMiddle, yN = N;
    break;
  case 3:
    x0 = xMiddle, xM = M, y0 = 0, yN = N - yMiddle;
    break;
  }
  std::array<int, 4> limits{x0, xM, y0, yN};
  return limits;
}

template <std::size_t M, std::size_t N>
void debugSendPrint(const std::array<std::array<int, M>, N> &grid,
                    const int currRank, const int nextRank, int x0, int xM,
                    int y0, int yN) {
  auto message =
      std::string(">> Process " + std::to_string(currRank) + " doing:\n");
  switch (currRank) {
  case 0:
    message += "  Top Left Sub-grid:\n";
    break;
  case 1:
    message += "  Bottom Left Sub-grid:\n";
    break;
  case 2:
    message += "  Bottom Right Sub-grid:\n";
    break;
  case 3:
    message += "  Top Right Sub-grid:\n";
    break;
  }
  for (int i = y0; i < yN; ++i) {
    for (int j = x0; j < xM; ++j) {
      message += std::to_string(grid[i][j]);
      message += "|";
    }
    message += '\n';
  }
  message += "  And sending it to process " + std::to_string(nextRank) + "\n";
  message += "----------\n";
  std::cout << message;
}

void debugReceivePrint(int storage[], const int currRank, const int prevRank,
                       std::pair<int, int> buffSize) {
  auto message = std::string("<< Process ");
  message = message + std::to_string(currRank) + " received from " +
            std::to_string(prevRank) + "\n";
  for (int i = 0; i < buffSize.first * buffSize.second; i++) {
    if (i == 0) {
      message = message + " ";
    } else if (i % buffSize.first == 0) {
      message = message + "\n ";
    }
    message = message + std::to_string(storage[i]) + "| ";
  }
  message = message + "\n==========\n";
  std::cout << message;
}

template <std::size_t W, std::size_t L>
void Send(std::array<std::array<int, W>, L> &grid, int rank, int nextRank,
          int M, int N) {
  auto [x0, xM, y0, yN] = GetLimits(rank, M, N);
  debugSendPrint(grid, rank, nextRank, x0, xM, y0, yN);
  auto tmpBuf = prepareSubGrid(grid, x0, xM, y0, yN);
  std::pair<int, int> buffSize(xM - x0, yN - y0);
  MPI_Send(&buffSize, 2, MPI_INT, nextRank, 1, MPI_COMM_WORLD);
  MPI_Send(&tmpBuf[0], buffSize.first * buffSize.second, MPI_INT, nextRank, 0,
           MPI_COMM_WORLD);
  delete[] tmpBuf;
}

void Receive(int rank, int prevRank, int M, int N) {
  std::pair<int, int> bS;
  MPI_Recv(&bS, 2, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  int storage[bS.first * bS.second];
  MPI_Recv(&storage, bS.first * bS.second, MPI_INT, prevRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  debugReceivePrint(storage, rank, prevRank, bS);
}

const int M = 5, N = 5;
std::array<std::array<int, M>, N> grid = {1, 1, 1, 1, 9, 2, 2, 2, 2, 8, 3, 3, 3,
                                          3, 7, 4, 4, 4, 4, 6, 5, 5, 5, 5, 0};

int main(int argc, char **argv) {
  int rank, size;
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Fill data for grid only in process 0
  if (rank == 3) {
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
    Send(grid, rank, nextRank, M, N);
    Receive(rank, prevRank, M, N);
  } else {
    Receive(rank, prevRank, M, N);
    Send(grid, rank, nextRank, M, N);
  }
  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}