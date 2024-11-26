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

void ExchangeResid(eDir dir, int nextRank, int prevRank) {

};
void ExchangeTau();
void ExchangeSolut();
void ExchangeMaxDiff();

// const int M = 40, N = 40;
// const int maxIter = 1e5;
const int M = 6, N = 6;
const int maxIter = 4;
const double tolerance = 1e-6;
const double h1 = 3.0 / M, h2 = 3.0 / N;
// auto method = sMethod::lin;
int masterRank = 0;
int tagData = 40, tagTau = 50, tagResid = 60, tagMaxDiff = 70, tagSol = 80;
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

  auto [x0, xM, y0, yN] = GetSectors(size, rank, M, N);
  auto domainSolution =
      Solution(xM - x0, yN - y0, x0, y0, h1, h2, maxIter, tolerance);
  domainSolution.ComputeABF();
  std::pair<int, int> sizeAndState(0, 0);
  std::pair<double, double> tauNomDenom(0, 0);
  double tau = 0.0;
  double maxDiff;

  if (size == 2) {

    for (int iter = 0; iter < maxIter; iter++) {
      maxDiff = 0.0;
      if (rank % 2 == 0) {
        // Calculate and exchange residuals
        domainSolution.CalculateResid();
        auto r = domainSolution.GetResidBoarder(eDir::right);
        MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, tagResid, comm);
        MPI_Recv(r.data(), r.size(), MPI_DOUBLE, prevRank, tagResid, comm,
                 status);
        domainSolution.SetResidBoarder(right, r);
        // Calculate and exchange steps tau
        auto tauPart = domainSolution.CalculateTau();
        MPI_Send(&tauPart, 2, MPI_DOUBLE, nextRank, tagTau, comm);
        MPI_Recv(&tauNomDenom, 2, MPI_DOUBLE, prevRank, tagTau, comm, status);
        auto [n1, d1] = tauPart;
        auto [n2, d2] = tauNomDenom;
        tau = (n1 + n2) / (d1 + d2);
      } else {
        // Calculate and exchange residuals
        domainSolution.CalculateResid();
        Grid::line_t r(domainSolution.GetN() + 1);
        MPI_Recv(r.data(), r.size(), MPI_DOUBLE, prevRank, tagResid, comm,
                 status);
        domainSolution.SetResidBoarder(left, r);
        r = domainSolution.GetResidBoarder(eDir::left);
        MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, tagResid, comm);
        // Calculate and exchange steps tau
        auto tauPart = domainSolution.CalculateTau();
        MPI_Recv(&tauNomDenom, 2, MPI_DOUBLE, prevRank, tagTau, comm, status);
        MPI_Send(&tauPart, 2, MPI_DOUBLE, nextRank, tagTau, comm);
        auto [n1, d1] = tauPart;
        auto [n2, d2] = tauNomDenom;
        tau = (n1 + n2) / (d1 + d2);
        // tau += tauPart;
      }
      // printf("rank %d, iter %d, tau=%f\n", rank, iter, tau);
      // At this point residuals are exchanged and step is common for both ranks
      MPI_Barrier(comm);

      if (rank % 2 == 0) {
        // Calculate and exchange solutions
        domainSolution.OneStepOfSolution(tau);
        auto w = domainSolution.GetSolutBoarder(right);
        MPI_Send(w.data(), w.size(), MPI_DOUBLE, nextRank, tagSol, comm);
        MPI_Recv(w.data(), w.size(), MPI_DOUBLE, prevRank, tagSol, comm,
                 status);
        domainSolution.SetSolutBoarder(right, w);
        // Calculate and exchange maximum of norm
        auto diff = domainSolution.CalculateMaxDiff(tau);
        MPI_Send(&diff, 1, MPI_DOUBLE, nextRank, tagMaxDiff, comm);
        MPI_Recv(&maxDiff, 1, MPI_DOUBLE, prevRank, tagMaxDiff, comm, status);
        maxDiff += diff;
      } else {
        domainSolution.OneStepOfSolution(tau);
        Grid::line_t w(domainSolution.GetN() + 1);
        MPI_Recv(w.data(), w.size(), MPI_DOUBLE, prevRank, tagSol, comm,
                 status);
        domainSolution.SetSolutBoarder(left, w);
        w = domainSolution.GetSolutBoarder(eDir::left);
        MPI_Send(w.data(), w.size(), MPI_DOUBLE, nextRank, tagSol, comm);
        // Calculate and exchange maximum of norm
        auto diff = domainSolution.CalculateMaxDiff(tau);
        MPI_Recv(&maxDiff, 1, MPI_DOUBLE, prevRank, tagMaxDiff, comm, status);
        MPI_Send(&diff, 1, MPI_DOUBLE, nextRank, tagMaxDiff, comm);
        maxDiff += diff;
      }
      printf("rank %d, iter %d, maxDiff=%f\n", rank, iter, maxDiff);
      MPI_Barrier(comm);

      if (maxDiff < tolerance || iter == maxIter - 1) {
        if (rank % 2 == 0) {
          break;
        } else {
          auto flattened = domainSolution.Flatten(left, 2);
          MPI_Send(flattened.data(), flattened.size(), MPI_DOUBLE, masterRank,
                   tagData, comm);
          break;
        }
      }
    }

    if (rank == masterRank) {
      // Grid::line_t boardVals = ReceiveMaster(sizeAndState, prevRank);
      auto size =
          ((domainSolution.GetM() + 1) - 2) * (domainSolution.GetN() + 1);
      Grid::line_t solVals(size);
      MPI_Recv(solVals.data(), solVals.size(), MPI_DOUBLE, prevRank, tagData,
               comm, status);
      //  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // domainSolution.Print();
      auto joinedGrid = domainSolution.Join(solVals, eDir::right, 2);
      // std::cout << "Flattened:\n";
      // for (auto elem: boardVals){
      //   std::cout << elem << ", ";
      // }
      // std::cout << std::endl;
      // std::cout << "\nSolution after receiving right boarder:\n" <<
      // joinedGrid.Print();
      joinedGrid.SaveToFile("mpi");
    }
  } else {

    for (int iter = 0; iter < maxIter; ++iter) {
      if (rank % 2 == 0) {
        // Calculate and exchange residuals
        domainSolution.CalculateResid();
        Grid::line_t r;
        Grid::line_t r2;
        // eDir dir = rank == 0 ? bottom : top;
        if (rank == 0) {
          r = domainSolution.GetResidBoarder(eDir::bottom);
          r2 = domainSolution.GetResidBoarder(eDir::right);
        } else {
          r = domainSolution.GetResidBoarder(eDir::top);
          r2 = domainSolution.GetResidBoarder(eDir::left);
        }
        MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, tagResid, comm);
        MPI_Send(r2.data(), r2.size(), MPI_DOUBLE, prevRank, tagResid, comm);
        MPI_Recv(r.data(), r.size(), MPI_DOUBLE, prevRank, tagResid, comm,
                 status);
        MPI_Recv(r2.data(), r2.size(), MPI_DOUBLE, nextRank, tagResid, comm,
                 status);
        if (rank == 0) {
          domainSolution.SetResidBoarder(right, r);
          domainSolution.SetResidBoarder(bottom, r2);
        } else {
          domainSolution.SetResidBoarder(left, r);
          domainSolution.SetResidBoarder(top, r2);
        }

      } else {
        // Calculate and exchange residuals
        domainSolution.CalculateResid();
        Grid::line_t r(domainSolution.GetN() + 1);
        Grid::line_t r2(domainSolution.GetN() + 1);
        MPI_Recv(r.data(), r.size(), MPI_DOUBLE, prevRank, tagResid, comm,
                 status);
        MPI_Recv(r2.data(), r2.size(), MPI_DOUBLE, nextRank, tagResid, comm,
                 status);
        if (rank == 1) {
          domainSolution.SetResidBoarder(top, r);
          r = domainSolution.GetResidBoarder(right);
          domainSolution.SetResidBoarder(right, r2);
          r2 = domainSolution.GetResidBoarder(top);
        } else {
          domainSolution.SetResidBoarder(bottom, r);
          r = domainSolution.GetResidBoarder(left);
          domainSolution.SetResidBoarder(left, r2);
          r2 = domainSolution.GetResidBoarder(bottom);
        }
        MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, tagResid, comm);
        MPI_Send(r2.data(), r2.size(), MPI_DOUBLE, prevRank, tagResid, comm);
      }
      // printf("rank %d, iter %d, tau=%f\n", rank, iter, tau);
      // At this point residuals are exchanged and step is common for both ranks
      MPI_Barrier(comm);
    }
    domainSolution.Print();
    if (rank == masterRank) {
    }
  }
  MPI_Finalize();
  return 0;
}