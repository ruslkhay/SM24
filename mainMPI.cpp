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
  std::cout << message << std::endl;
}

const int M = 5, N = 5;
std::array<std::array<int, M>, N> grid = {1, 1, 1, 1, 9, 2, 2, 2, 2, 8, 3, 3, 3,
                                          3, 7, 4, 4, 4, 4, 6, 5, 5, 5, 5, 0};
// for (int i = 0; i < M; i++) {
//   for (int j = 0; j < N; j++) {
//     grid[i][j] = rand() % std::max(M, N); // Random values
//   }
// }

int main(int argc, char **argv) {
  int rank, size;
  const int xMiddle = M / 2 + 1, yMiddle = N / 2 + 1;
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

  // Distribute the data accordingly
  if (rank == 0) {
    int nextRank = 1;
    int x0 = 0, xM = xMiddle, y0 = 0, yN = N - yMiddle;
    debugSendPrint(grid, rank, nextRank, x0, xM, y0, yN);
    auto tmpBuf = prepareSubGrid(grid, x0, xM, y0, yN);
    int count = (xM - x0) * (yN - y0);
    MPI_Send(&tmpBuf[0], count, MPI_INT, nextRank, 0, MPI_COMM_WORLD);
    delete[] tmpBuf;
  } else if (rank == 1) {
    int nextRank = 0;
    int x0 = 0, xM = xMiddle, y0 = N - yMiddle, yN = N;
    debugSendPrint(grid, rank, nextRank, x0, xM, y0, yN);
    auto tmpBuf = prepareSubGrid(grid, x0, xM, y0, yN);
    int count = (xM - x0) * (yN - y0);
    MPI_Send(&tmpBuf[0], count, MPI_INT, nextRank, 0, MPI_COMM_WORLD);
    delete[] tmpBuf;
  } else if (rank == 2) {
    int nextRank = 2;
    int x0 = xMiddle, xM = M, y0 = N - yMiddle, yN = N;
    debugSendPrint(grid, rank, nextRank, x0, xM, y0, yN);
    auto tmpBuf = prepareSubGrid(grid, x0, xM, y0, yN);
    int count = (xM - x0) * (yN - y0);
    MPI_Send(&tmpBuf[0], count, MPI_INT, nextRank, 0, MPI_COMM_WORLD);
    delete[] tmpBuf;
  } else if (rank == 3) {
    int nextRank = 3;
    int x0 = xMiddle, xM = M, y0 = 0, yN = N - yMiddle;
    debugSendPrint(grid, rank, nextRank, xMiddle, M, 0, N - yMiddle);
    auto tmpBuf = prepareSubGrid(grid, xMiddle, M, 0, N - yMiddle);
    int count = (xM - x0) * (yN - y0);
    MPI_Send(&tmpBuf[0], count, MPI_INT, nextRank, 0, MPI_COMM_WORLD);
    delete[] tmpBuf;
  }

  // Now each process should receive its respective sub-grid
  if (rank == 1) {
    int prevRank = 0;
    int first[2][3]; // Process 0
    MPI_Recv(&first, 2 * 3, MPI_INT, prevRank, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("<< Process %d received:\n", rank);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%4d", first[i][j]);
      }
      printf("\n");
    }
    std::cout << "----------\n";
  } else if (rank == 0) {
    int prevRank = 1;
    int second[3][3]; // Process 1
    MPI_Recv(&second, 3 * 3, MPI_INT, prevRank, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("<< Process %d received:\n", rank);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        printf("%4d", second[i][j]);
      }
      printf("\n");
    }
    std::cout << "----------\n";
  } else if (rank == 2) {
    int prevRank = 2;
    int third[3][2];
    MPI_Recv(&third, 3 * 2, MPI_INT, prevRank, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("<< Process %d received:\n", rank);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        printf("%4d", third[i][j]);
      }
      printf("\n");
    }
    std::cout << "----------\n";
  } else if (rank == 3) {
    int prevRank = 3;
    int fourth[2][2];
    MPI_Recv(&fourth, 2 * 2, MPI_INT, prevRank, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("<< Process %d received:\n", rank);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        printf("%4d", fourth[i][j]);
      }
      printf("\n");
    }
    std::cout << "----------\n";
  }
  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}