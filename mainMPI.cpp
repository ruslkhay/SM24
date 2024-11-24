#include "src/mpi.hpp"
#include "src/solution.hpp"
#include <algorithm>
#include <cmath>
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

void _Receive(int rank, int prevRank) {
  std::pair<int, int> bS;
  MPI_Recv(&bS, 2, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::vector<double> storage(bS.first);
  MPI_Recv(&storage[0], bS.first, MPI_DOUBLE, prevRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  // debugReceivePrint(storage, rank, prevRank, bS);
}

Grid::line_t Receive(int rank, int prevRank) {
  std::pair<int, int> bS;
  MPI_Recv(&bS, 2, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  Grid::line_t storage(bS.first);
  MPI_Recv(&storage[0], bS.first, MPI_DOUBLE, prevRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  // debugReceivePrint(storage, rank, prevRank, bS);
  return storage;
}

const int M = 40, N = 40;
// const int maxIter = 1e5;
// const int M = 4, N = 4;
const int maxIter = 801;
const double tolerance = 1e-6;
const double h1 = 3.0 / M, h2 = 3.0 / N;
auto method = sMethod::lin;

void SendToMaster(const std::pair<int, int> sizeAndState,
                  const Grid::line_t &boardVals, int masterRank,
                  int tagSAS = 665, int tagData = 666) {
  MPI_Send(&sizeAndState, 2, MPI_INT, masterRank, tagSAS, MPI_COMM_WORLD);
  MPI_Send(&boardVals[0], sizeAndState.first, MPI_DOUBLE, masterRank, tagData,
           MPI_COMM_WORLD);
}
Grid::line_t ReceiveByMaster(std::pair<int, int> &sizeAndState, int sender,
                             int tagSAS = 665, int tagData = 666) {
  MPI_Recv(&sizeAndState, 2, MPI_INT, sender, tagSAS, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  Grid::line_t boardVals(sizeAndState.first);
  MPI_Recv(&boardVals[0], sizeAndState.first, MPI_DOUBLE, sender, tagData,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  return boardVals;
}

void Send(const std::pair<int, int> sizeAndState, const Grid::line_t &boardVals,
          const std::pair<double, double> tauNomDenom, int nextRank,
          int tagSAS = 0, int tagData = 1, int tagTau = 2) {
  MPI_Send(&sizeAndState, 2, MPI_INT, nextRank, tagSAS, MPI_COMM_WORLD);
  MPI_Send(&boardVals[0], sizeAndState.first, MPI_DOUBLE, nextRank, tagData,
           MPI_COMM_WORLD);
  MPI_Send(&tauNomDenom, 2, MPI_DOUBLE, nextRank, tagTau, MPI_COMM_WORLD);
}

Grid::line_t Receive(std::pair<int, int> &sizeAndState,
                     std::pair<double, double> &tauNomDenom, int prevRank,
                     int tagSAS = 0, int tagData = 1, int tagTau = 2) {
  MPI_Recv(&sizeAndState, 2, MPI_INT, prevRank, tagSAS, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  Grid::line_t boardVals(sizeAndState.first);
  MPI_Recv(&boardVals[0], sizeAndState.first, MPI_DOUBLE, prevRank, tagData,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Get tau
  MPI_Recv(&tauNomDenom, 2, MPI_DOUBLE, prevRank, tagTau, MPI_COMM_WORLD,
           MPI_STATUSES_IGNORE);
  return boardVals;
}

int main(int argc, char **argv) {
  int rank, size;
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::pair<int, int> buffSize(0, 0);
  int nextRank = (rank + 1) % size;
  int prevRank = rank == 0 ? size - 1 : rank - 1;
  int masterRank = 0;

  auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
  auto domainSolution =
      Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
  domainSolution.ComputeABF();
  std::pair<int, int> sizeAndState(0, 0);

  std::pair<double, double> tauNomDenom(0, 0);
  double tau = 0.0;

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank % 2 == 0) {
    for (int iter = 0; iter < maxIter; iter++) {
      // Store maximum of norm ||w_(k+1) - w_k||
      auto maxDiff = 0.0;
      // Produce one step and store ||w_(k+1) - w_k|| into diff

      // if (iter == 0 || sizeAndState.second) {
      if (iter == 0) {
        tauNomDenom = domainSolution.CalculateTau();
        tau = tauNomDenom.first / tauNomDenom.second;
      } else {
        auto [tNom, tDenom] = domainSolution.CalculateTau();
        tau = (tauNomDenom.first + tNom) / (tauNomDenom.second + tDenom);
      }

      auto diff = domainSolution.OneStepOfSolution(tau);
      maxDiff = std::max(maxDiff, diff);
      // Send boarder values to next process
      // Check if builded solution is suitable for current domain
      auto rightBoardVals = domainSolution.GetColumn(domainSolution.GetM() - 1);

      // debugSendPrint(domainSolution.GetNodes(), rank, nextRank, x0, xM, y0,
      // yN,
      //                rightBoardVals);
      printf("maxDiff=%.8f, tau=%f for rank %d, iter№ %d, prev finished %d\n",
             maxDiff, tau, rank, iter, sizeAndState.second);

      if ((maxDiff < tolerance) || iter == maxIter - 1) {
        if (!sizeAndState.second) {
          sizeAndState = {rightBoardVals.size(), 1};
          Send(sizeAndState, rightBoardVals, tauNomDenom, nextRank);
        }
        break; // it is masterRank
        // sizeAndState = {rightBoardVals.size(), 1};
        // Send(sizeAndState, rightBoardVals, tauNomDenom, nextRank);
      }
      if (!sizeAndState.second) {
        // Send size, state (finish or not), and boarder values
        sizeAndState = {rightBoardVals.size(), 0};
        Send(sizeAndState, rightBoardVals, tauNomDenom, nextRank);
        auto boardVals = Receive(sizeAndState, tauNomDenom, prevRank);
        domainSolution.SetRightBoarder(boardVals);
      }
    }
  } else {
    for (int iter = 0; iter < maxIter; iter++) {
      // Store maximum of norm ||w_(k+1) - w_k||
      auto maxDiff = 0.0;
      // if process 0 is not finished
      if (!sizeAndState.second) { // for initial iteration it is true
        auto boardVals = Receive(sizeAndState, tauNomDenom, prevRank);
        // Add boarder values to domain
        domainSolution.SetLeftBoarder(boardVals);

        auto [tNom, tDenom] = domainSolution.CalculateTau();
        tau = (tauNomDenom.first + tNom) / (tauNomDenom.second + tDenom);
        tauNomDenom = {tNom, tDenom}; // This is for sending it to proc 0
      } else {
        tauNomDenom = domainSolution.CalculateTau();
        tau = tauNomDenom.first / tauNomDenom.second;
      }

      auto diff = domainSolution.OneStepOfSolution(tau);
      maxDiff = std::max(maxDiff, diff);
      printf("maxDiff=%.8f, tau=%f for rank %d, iter№ %d, prev finished %d\n",
             maxDiff, tau, rank, iter, sizeAndState.second);
      // Check if builded solution is suitable for current domain
      if (!sizeAndState.second) {
        auto leftBoardVals = domainSolution.GetColumn(1);
        sizeAndState = {leftBoardVals.size(), 0};
        // debugSendPrint(domainSolution.GetNodes(), rank, nextRank, x0, xM, y0,
        //                yN, leftBoardVals);
        Send(sizeAndState, leftBoardVals, tauNomDenom, nextRank);
      }
      if (maxDiff < tolerance || iter == maxIter - 1) {
        // std::cout << "Before flattening:\n";
        if (!sizeAndState.second) {
          auto leftBoardVals = domainSolution.GetColumn(1);
          sizeAndState = {leftBoardVals.size(), 1};
          // debugSendPrint(domainSolution.GetNodes(), rank, nextRank, x0, xM,
          // y0,
          //                yN, leftBoardVals);
          Send(sizeAndState, leftBoardVals, tauNomDenom, nextRank);
        }
        auto flattened = domainSolution.Flatten(eDir::left, 2);
        sizeAndState = {flattened.size(), 1};
        SendToMaster(sizeAndState, flattened, masterRank);
        // Send(sizeAndState, flattened, tauNomDenom, nextRank);
        break;
      }
    }
  }

  if (rank == masterRank) {
    auto boardVals = ReceiveByMaster(sizeAndState, prevRank);
    // domainSolution.Print();
    auto joinedGrid = domainSolution.Join(boardVals, eDir::right, 2);
    // std::cout << "Flattened:\n";
    // for (auto elem: boardVals){
    //   std::cout << elem << ", ";
    // }
    // std::cout << std::endl;
    // std::cout << "Solution after receiving right boarder:\n" << std::endl;
    // joinedGrid.Print();
    joinedGrid.SaveToFile("mpi");
  }

  MPI_Finalize();
  return 0;
}