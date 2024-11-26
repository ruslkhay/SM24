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

int tagSAS = 0, tagData = 1, tagTau = 3;
void Send(const std::pair<int, int> sizeAndState,
          const std::pair<Grid::line_t, Grid::line_t> &boardVals,
          const std::pair<double, double> tauNomDenom, int nextRank) {
  // auto size = sizeAndState.first;
  MPI_Send(&sizeAndState, 2, MPI_INT, nextRank, tagSAS, MPI_COMM_WORLD);
  // Send solution and residuals
  auto w = boardVals.first;
  MPI_Send(&w[0], w.size(), MPI_DOUBLE, nextRank, tagData, MPI_COMM_WORLD);
  auto r = boardVals.second;
  MPI_Send(&r[0], r.size(), MPI_DOUBLE, nextRank, tagData + 1, MPI_COMM_WORLD);
  // Send step tau nominator and denominator
  MPI_Send(&tauNomDenom, 2, MPI_DOUBLE, nextRank, tagTau, MPI_COMM_WORLD);
}

void SendMaster(const std::pair<int, int> sizeAndState,
                const Grid::line_t &boardVals, int nextRank) {
  MPI_Send(&sizeAndState, 2, MPI_INT, nextRank, tagSAS, MPI_COMM_WORLD);
  // Send solution and residuals
  MPI_Send(&boardVals[0], sizeAndState.first, MPI_DOUBLE, nextRank, tagData,
           MPI_COMM_WORLD);
}

Grid::line_t ReceiveMaster(std::pair<int, int> &sizeAndState, int prevRank) {
  MPI_Recv(&sizeAndState, 2, MPI_INT, prevRank, tagSAS, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  Grid::line_t boardVals(sizeAndState.first);
  MPI_Recv(&boardVals[0], sizeAndState.first, MPI_DOUBLE, prevRank, tagData,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  return boardVals;
}

std::pair<Grid::line_t, Grid::line_t>
Receive(std::pair<int, int> &sizeAndState,
        std::pair<double, double> &tauNomDenom, int prevRank) {
  MPI_Recv(&sizeAndState, 2, MPI_INT, prevRank, tagSAS, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  // Solution receive
  Grid::line_t w(sizeAndState.first);
  MPI_Recv(&w[0], sizeAndState.first, MPI_DOUBLE, prevRank, tagData,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Residuals receive
  Grid::line_t r(sizeAndState.first);
  MPI_Recv(&r[0], sizeAndState.first, MPI_DOUBLE, prevRank, tagData + 1,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Step tau
  MPI_Recv(&tauNomDenom, 2, MPI_DOUBLE, prevRank, tagTau, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  // std::pair<Grid::line_t, Grid::line_t> boardVals;
  // return boardVals;
  return {w, r};
}

const int M = 4, N = 4;
// const int maxIter = 1e5;
// const int M = 6, N = 6;
const int maxIter = 8;
const double tolerance = 1e-6;
const double h1 = 3.0 / M, h2 = 3.0 / N;
// auto method = sMethod::lin;

int tagResid = 6;
auto comm = MPI_COMM_WORLD;
auto status = MPI_STATUS_IGNORE;
int main(int argc, char **argv) {
  int rank, size;
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::pair<int, int> buffSize(0, 0);
  int nextRank = (rank + 1) % size;
  int prevRank = rank == 0 ? size - 1 : rank - 1;
  // int masterRank = 0;

  auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
  auto domainSolution =
      Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
  domainSolution.ComputeABF();
  std::pair<int, int> sizeAndState(0, 0);
  std::pair<double, double> tauNomDenom(0, 0);
  double tau = 0.0;

  for (int iter = 0; iter < maxIter; iter++) {

    if (rank % 2 == 0) {
      // Calculate and exchange residuals;
      domainSolution.CalculateResid();
      auto r = domainSolution.GetResidBoarder(eDir::right);
      MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, 6, comm);
      MPI_Recv(r.data(), r.size(), MPI_DOUBLE, prevRank, 6, comm, status);
      domainSolution.SetResidBoarder(right, r);
      // Calculate and exchange steps tau
      auto [tNom, tDenom] = domainSolution.CalculateTau();
      auto tauPart = tNom / tDenom;
      MPI_Send(&tauPart, 1, MPI_DOUBLE, nextRank, tagTau, comm);
      MPI_Recv(&tau, 1, MPI_DOUBLE, prevRank, tagTau, comm, status);
      tau += tauPart;
    } else {
      domainSolution.CalculateResid();
      Grid::line_t r(domainSolution.GetN() + 1);
      MPI_Recv(r.data(), r.size(), MPI_DOUBLE, prevRank, 6, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      domainSolution.SetResidBoarder(left, r);
      r = domainSolution.GetResidBoarder(eDir::left);
      MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, 6, MPI_COMM_WORLD);
      // Calculate and exchange steps tau
      auto [tNom, tDenom] = domainSolution.CalculateTau();
      auto tauPart = tNom / tDenom;
      MPI_Recv(&tau, 1, MPI_DOUBLE, prevRank, tagTau, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      tau += tauPart;
      MPI_Send(&tauPart, 1, MPI_DOUBLE, nextRank, tagTau, MPI_COMM_WORLD);
    }
    printf("Rank №%d tau = %f\n", rank, tau);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // if (rank == masterRank) {
  //   Grid::line_t boardVals = ReceiveMaster(sizeAndState, prevRank);
  //   //  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   // domainSolution.Print();
  //   auto joinedGrid = domainSolution.Join(boardVals, eDir::right, 2);
  //   // std::cout << "Flattened:\n";
  //   // for (auto elem: boardVals){
  //   //   std::cout << elem << ", ";
  //   // }
  //   // std::cout << std::endl;
  //   // std::cout << "\nSolution after receiving right boarder:\n" <<
  //   std::endl;
  //   // joinedGrid.Print();
  //   // joinedGrid.SaveToFile("mpi");
  // }

  MPI_Finalize();
  return 0;
}