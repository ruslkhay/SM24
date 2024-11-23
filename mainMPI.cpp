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

std::vector<std::vector<double>>
SumUpCols(const std::vector<std::vector<double>> &matrix, int firstCol) {
  // Initialize the resulting matrix
  std::vector<std::vector<double>> result;

  // Check if the matrix is empty or firstCol is out of bounds
  if (matrix.empty() || static_cast<int>(matrix[0].size()) <= firstCol) {
    return result; // Return an empty matrix in this case
  }

  // Iterate through each row of the matrix
  for (const auto &row : matrix) {
    std::vector<double> newRow;

    // Add elements from the original row except those at firstCol + 1
    for (int i = 0; i < static_cast<int>(row.size()); ++i) {
      if (i == firstCol) {
        newRow.push_back(row[i]); // Keep the first column as it is
      } else if (i == firstCol + 1) {
        // Sum up the values of the two columns
        newRow.push_back(row[firstCol] + row[i]);
      } else {
        // Add elements from other columns
        newRow.push_back(row[i]);
      }
    }

    result.push_back(newRow);
  }

  // Resize the result to keep only the needed columns
  for (auto &row : result) {
    if (static_cast<int>(row.size()) > firstCol + 1) {
      row.erase(row.begin() + firstCol + 1); // Remove column at firstCol + 1
    }
  }

  return result;
}

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

const int M = 4, N = 4;
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
    // Get limits of first domain (0)
    auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
    auto solution =
        Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
    solution.Find(method);
    // auto rightGrid = solution.Flatten(eDir::left);

    // auto boardVals = solution.GetRightBoarder();
    // Take right inner nodes and send it to next domain
    auto boardVals = solution.GetColumn(solution.GetM() - 1);
    debugSendPrint(solution.GetNodes(), rank, nextRank, x0, xM, y0, yN,
                   boardVals);
    Send(boardVals, size, rank, nextRank);
    // Receive all nodes from next domain and join to solution
    auto boarderVal = Receive(rank, prevRank);
    auto joinedGrid = solution.Join(boarderVal, eDir::right);
    std::cout << "Solution after receiving right boarder:\n" << std::endl;
    joinedGrid.Print();

  } else {
    auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
    printf("(%d; %d), (%d, %d)\n", x0, xM, y0, yN);
    auto solution =
        Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
    auto boarderVal = Receive(rank, prevRank);
    solution.SetLeftBoarder(boarderVal);
    solution.Find(method);
    auto flattened = solution.Flatten(eDir::left);
    // auto boardVals = solution.GetRightBoarder();
    // auto boardVals = solution.GetColumn(solution.GetM() - 1);
    // debugSendPrint(solution.GetNodes(), rank, nextRank, x0, xM, y0, yN,
    //                boardVals);
    debugSendPrint(solution.GetNodes(), rank, nextRank, x0, xM, y0, yN,
                   flattened);
    Send(flattened, size, rank, nextRank);
  }
  // Finalize the MPI environment
  MPI_Finalize();
  return 0;
}