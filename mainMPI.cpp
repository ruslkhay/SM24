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
// // function to initialize and returning array
// template<std::size_t R1, std::size_t C1>
// void subGrid(std::array<std::array<int, R1>, C1>& grid ,int x0, int xM, int
// y0, int yN){
//     int** arr = new int*[yN-y0];
//     int N = yN - y0, M = xM - x0;
//     for (int i = x0; i < xM; ++i) {
//         arr[i] = new int[yN-y0];
//         for (int j = y0; j < yN; ++j) {
//             arr[i][j] = grid[i][j];
//         }
//     }
//     return arr;
// }

template <std::size_t M, std::size_t N>
int *prepareSubGrid(std::array<std::array<int, M>, N> &grid, int x0, int xM,
                    int y0, int yN) {
  int count = (xM - x0) * (yN - y0);
  int *tmpBuf = new int[count];
  for (int i = x0; i < xM; ++i) {
    for (int j = y0; j < yN; ++j) {
      tmpBuf[i * yN + j] = grid[i][j];
    }
  }
  return tmpBuf;
}

const int M = 5, N = 5;
std::array<std::array<int, M>, N> grid = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3,
                                          3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5};
// for (int i = 0; i < M; i++) {
//   for (int j = 0; j < N; j++) {
//     grid[i][j] = rand() % std::max(M, N); // Random values
//   }
// }

int main(int argc, char **argv) {
  int rank, size;

  // Declare the grid that will only be filled by rank 0
  // int grid[M][N];

  const int topLeftRows = 2, topLeftCols = 3;
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
  int first[2][3];  // Process 0
  int second[3][3]; // Process 1
  // int third[2][2];  // Process 2
  // int fourth[3][2]; // Process 3

  // Distribute the data accordingly
  if (rank == 0) {
    // printf("On process %d; processing ", rank);
    // Print Top Left Sub-grid
    std::cout << "Top Left Sub-grid:\n";
    for (int i = 0; i < topLeftRows; ++i) {
      for (int j = 0; j < topLeftCols; ++j) {
        std::cout << grid[i][j] << "|";
      }
      std::cout << "\n";
    }
    // Process 0: top-left 2x3 grid
    int count = topLeftRows * topLeftCols;
    int tmpBuf[count];
    for (int i = 0; i < topLeftRows; ++i) {
      for (int j = 0; j < topLeftCols; ++j) {
        tmpBuf[i * topLeftCols + j] = grid[i][j];
      }
    }
    MPI_Send(&tmpBuf[0], count, MPI_INT, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 1) {
    // printf("On process %d; processing ", rank);
    // Print Bottom Left Sub-grid
    std::cout << "Bottom Left Sub-grid:\n";
    for (int i = topLeftRows; i < topLeftRows + bottomLeftRows; ++i) {
      for (int j = 0; j < bottomLeftCols; ++j) {
        std::cout << grid[i][j] << "|";
      }
      std::cout << "\n";
    }
    // Process 1: bottom-left 3x3 grid
    // auto tmpBuf = prepareSubGrid(grid, 0, topRightRows, 0, topRightCols);
    int count = topRightRows * topRightCols;
    int tmpBuf[count];
    for (int i = 0; i < topRightRows; ++i) {
      for (int j = 0; j < topRightCols; ++j) {
        tmpBuf[i * topRightCols + j] = grid[i][j];
      }
    }
    MPI_Send(&tmpBuf[0], 3 * 3, MPI_INT, 1, 0, MPI_COMM_WORLD);
    // delete[] tmpBuf;
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
  if (rank == 0) {
    MPI_Recv(&first, 2 * 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // printf("Process %d received:\n", rank);
    // for (int i = 0; i < 2; i++) {
    //   for (int j = 0; j < 3; j++) {
    //     printf("%4d", first[i][j]);
    //   }
    //   printf("\n");
    // }
  } else if (rank == 1) {
    MPI_Recv(&second, 3 * 3, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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