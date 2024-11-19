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
      int a = grid[i][j];
      tmpBuf[(i - y0) * (xM - x0) + j] = a;
    }
  }
  return tmpBuf;
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

  const int topLeftRows = 2, topLeftCols = 3;
  const int xMiddle = M / 2 + 1, yMiddle = N / 2 + 1;
  const int topRightRows = 2, topRightCols = 2;
  const int bottomLeftRows = 3, bottomLeftCols = 3;
  const int bottomRightRows = 3; //, bottomRightCols = 2;

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

  // Update the grid with data specific for each rank
  int sub_grid[3][3]; // Maximum size for receiving

  // Create buffers for each process based on sub-grid sizes
  // int third[2][2];  // Process 2
  // int fourth[3][2]; // Process 3

  // Distribute the data accordingly
  if (rank == 0) {
    int nextRank = 1;
    printf(">> Process %d doing:\n", rank);
    // Print Top Left Sub-grid
    std::cout << "  Top Left Sub-grid:\n";
    for (int i = 0; i < N - yMiddle; ++i) {
      for (int j = 0; j < xMiddle; ++j) {
        std::cout << grid[i][j] << "|";
      }
      std::cout << "\n";
    }
    printf("  And sending it to process %d\n", nextRank);
    // Process 0: top-left 2x3 grid
    auto tmpBuf = prepareSubGrid(grid, 0, xMiddle, 0, N - yMiddle);
    int count = xMiddle * (N - yMiddle);
    MPI_Send(&tmpBuf[0], count, MPI_INT, nextRank, 0, MPI_COMM_WORLD);
    delete[] tmpBuf;
  } else if (rank == 1) {
    int nextRank = 0;
    printf(">> Process %d doing:\n", rank);
    // Print Bottom Left Sub-grid
    std::cout << "Bottom Left Sub-grid:\n";
    for (int i = topLeftRows; i < topLeftRows + bottomLeftRows; ++i) {
      for (int j = 0; j < bottomLeftCols; ++j) {
        std::cout << grid[i][j] << "|";
      }
      std::cout << "\n";
    }
    // Process 1: bottom-left 3x3 grid
    auto tmpBuf = prepareSubGrid(grid, 0, xMiddle, N - yMiddle, N);
    int count = xMiddle * yMiddle;
    MPI_Send(&tmpBuf[0], count, MPI_INT, nextRank, 0, MPI_COMM_WORLD);
    delete[] tmpBuf;
  } else if (rank == 2) {
    // printf("On process %d; processing ", rank);
    // Print Bottom Right Sub-grid
    std::cout << "Bottom Right Sub-grid:\n";
    for (int i = topLeftRows; i < topLeftRows + bottomRightRows; ++i) {
      for (int j = topLeftCols; j < topLeftCols + topRightCols; ++j) {
        std::cout << grid[i][j] << "|";
      }
      std::cout << "\n";
    }
    // Process 2: top-right 2x2 grid
    MPI_Send(&grid[0][3], 2 * 2, MPI_INT, 2, 0, MPI_COMM_WORLD);
  } else if (rank == 3) {
    // printf("On process %d; processing ", rank);
    // Print Top Right Sub-grid
    std::cout << "Top Right Sub-grid:\n";
    for (int i = 0; i < topRightRows; ++i) {
      for (int j = topLeftCols; j < topLeftCols + topRightCols; ++j) {
        std::cout << grid[i][j] << "|";
      }
      std::cout << "\n";
    }
    // Process 3: bottom-right 3x2 grid
    MPI_Send(&grid[3][0], 3 * 2, MPI_INT, 3, 0, MPI_COMM_WORLD);
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
    // printf("Process %d received:\n", rank);
    // for (int i = 0; i < 3; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     printf("%4d", second[i][j]);
    //   }
    //   printf("\n");
    // }
  } else if (rank == 2) {
    MPI_Recv(&sub_grid, 2 * 2, MPI_INT, 2, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  } else if (rank == 3) {
    MPI_Recv(&sub_grid, 3 * 2, MPI_INT, 3, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }
  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}