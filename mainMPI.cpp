#include "src/solution.hpp"
#include <algorithm>
#include <chrono>
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

std::array<int, 4> GetLimitsTwoProc(int rank, int M, int N);
std::array<int, 4> GetLimitsFourProc(int rank, int M, int N);
std::array<int, 4> GetSectors(int procNum, int rank, int M, int N);

const int M = 40, N = 40;
const int maxIter = 1e5;
// const int M = 8, N = 8;
const double tolerance = 1e-6;
const double h1 = 3.0 / M, h2 = 3.0 / N;
// auto method = sMethod::lin;
int masterRank = 0;
int tagData = 40, tagTau = 50, tagResid = 60, tagMaxDiff = 70, tagSol = 80;
auto comm = MPI_COMM_WORLD;
auto status = MPI_STATUS_IGNORE;

int main(int argc, char **argv) {
  int rank, size;
  double start, stop;
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
  start = MPI_Wtime();
  domainSolution.ComputeABF();
  double tau = 0.0;
  double tauNom, tauDenom;
  double maxDiff;

  if (size == 1) {
    domainSolution.Find(lin);
    stop = MPI_Wtime();
    domainSolution._execTime = std::chrono::duration_cast<Solution::time_t>(
        std::chrono::duration<double>(stop - start));
    domainSolution.SaveToFile("mpi_1proc");
  } else if (size == 2) {
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
        auto [n, d] = tauPart;
        MPI_Allreduce(&n, &tauNom, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&d, &tauDenom, 1, MPI_DOUBLE, MPI_SUM, comm);
        tau = tauNom / tauDenom;

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
        auto [n, d] = domainSolution.CalculateTau();
        MPI_Allreduce(&n, &tauNom, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&d, &tauDenom, 1, MPI_DOUBLE, MPI_SUM, comm);
        tau = tauNom / tauDenom;
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
        MPI_Allreduce(&diff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, comm);
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
        MPI_Allreduce(&diff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, comm);
      }
      // printf("rank %d, iter %d, maxDiff=%f\n", rank, iter, maxDiff);
      // MPI_Barrier(comm);

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
      auto size =
          ((domainSolution.GetM() + 1) - 2) * (domainSolution.GetN() + 1);
      Grid::line_t solVals(size);
      MPI_Recv(solVals.data(), solVals.size(), MPI_DOUBLE, prevRank, tagData,
               comm, status);
      auto joinedGrid = domainSolution.Join(solVals, eDir::right, 2);
      stop = MPI_Wtime();
      joinedGrid._execTime = std::chrono::duration_cast<Solution::time_t>(
          std::chrono::duration<double>(stop - start));
      joinedGrid.SaveToFile("mpi_2proc");
    }
  } else if (size == 4) {
    //------------------------------------------------------------------------------
    for (int iter = 0; iter < maxIter; ++iter) {
      maxDiff = 0.0;
      if (rank % 2 == 0) {
        // Calculate and exchange solutions
        domainSolution.CalculateResid();
        Grid::line_t r;
        Grid::line_t r2;
        if (rank == 0) {
          r = domainSolution.GetResidBoarder(eDir::bottom);
          r2 = domainSolution.GetResidBoarder(eDir::right);
        } else {
          r = domainSolution.GetResidBoarder(eDir::top);
          r2 = domainSolution.GetResidBoarder(eDir::left);
        }
        r.resize(domainSolution.GetN() + 1);
        r2.resize(domainSolution.GetN() + 1);
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
        // Calculate and exchange tau
        auto [n, d] = domainSolution.CalculateTau();
        MPI_Allreduce(&n, &tauNom, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&d, &tauDenom, 1, MPI_DOUBLE, MPI_SUM, comm);
        tau = tauNom / tauDenom;
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
        r.resize(domainSolution.GetN() + 1);
        r2.resize(domainSolution.GetN() + 1);
        MPI_Send(r.data(), r.size(), MPI_DOUBLE, nextRank, tagResid, comm);
        MPI_Send(r2.data(), r2.size(), MPI_DOUBLE, prevRank, tagResid, comm);
        // Calculate and exchange steps tau
        auto [n, d] = domainSolution.CalculateTau();
        MPI_Allreduce(&n, &tauNom, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&d, &tauDenom, 1, MPI_DOUBLE, MPI_SUM, comm);
        tau = tauNom / tauDenom;
      }
      // printf("rank %d, iter %d, tau=%f\n", rank, iter, tau);
      // At this point residuals are exchanged and step is common for both ranks
      // MPI_Barrier(comm);
      if (rank % 2 == 0) {
        // Calculate and exchange residuals
        domainSolution.OneStepOfSolution(tau);
        Grid::line_t w;
        Grid::line_t w2;
        if (rank == 0) {
          w = domainSolution.GetSolutBoarder(eDir::bottom);
          w2 = domainSolution.GetSolutBoarder(eDir::right);
        } else {
          w = domainSolution.GetSolutBoarder(eDir::top);
          w2 = domainSolution.GetSolutBoarder(eDir::left);
        }
        w.resize(domainSolution.GetN() + 1);
        w2.resize(domainSolution.GetN() + 1);
        MPI_Send(w.data(), w.size(), MPI_DOUBLE, nextRank, tagSol, comm);
        MPI_Send(w2.data(), w2.size(), MPI_DOUBLE, prevRank, tagSol, comm);
        MPI_Recv(w.data(), w.size(), MPI_DOUBLE, prevRank, tagSol, comm,
                 status);
        MPI_Recv(w2.data(), w2.size(), MPI_DOUBLE, nextRank, tagSol, comm,
                 status);
        if (rank == 0) {
          domainSolution.SetSolutBoarder(right, w);
          domainSolution.SetSolutBoarder(bottom, w2);
        } else {
          domainSolution.SetSolutBoarder(left, w);
          domainSolution.SetSolutBoarder(top, w2);
        }
        // Calculate and exchange maximum of norm
        auto diff = domainSolution.CalculateMaxDiff(tau);
        MPI_Allreduce(&diff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, comm);
      } else {
        domainSolution.OneStepOfSolution(tau);
        Grid::line_t w(domainSolution.GetN() + 1);
        Grid::line_t w2(domainSolution.GetN() + 1);
        MPI_Recv(w.data(), w.size(), MPI_DOUBLE, prevRank, tagSol, comm,
                 status);
        MPI_Recv(w2.data(), w2.size(), MPI_DOUBLE, nextRank, tagSol, comm,
                 status);
        if (rank == 1) {
          domainSolution.SetSolutBoarder(top, w);
          w = domainSolution.GetSolutBoarder(right);
          domainSolution.SetSolutBoarder(right, w2);
          w2 = domainSolution.GetSolutBoarder(top);
        } else {
          domainSolution.SetSolutBoarder(bottom, w);
          w = domainSolution.GetSolutBoarder(left);
          domainSolution.SetSolutBoarder(left, w2);
          w2 = domainSolution.GetSolutBoarder(bottom);
        }
        w.resize(domainSolution.GetN() + 1);
        w2.resize(domainSolution.GetN() + 1);
        MPI_Send(w.data(), w.size(), MPI_DOUBLE, nextRank, tagSol, comm);
        MPI_Send(w2.data(), w2.size(), MPI_DOUBLE, prevRank, tagSol, comm);
        // Calculate and exchange steps tau
        auto diff = domainSolution.CalculateMaxDiff(tau);
        MPI_Allreduce(&diff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, comm);
      }
      // MPI_Barrier(comm);
      // printf("rank %d, iter %d, tau=%f, maxDiff=%f\n", rank, iter, tau,
      //  maxDiff);

      if (maxDiff < tolerance || iter == maxIter - 1) {
        if (rank == masterRank) {
          break;
        } else if (rank == 1) {
          auto flattened = domainSolution.Flatten(top, 2);
          MPI_Send(flattened.data(), flattened.size(), MPI_DOUBLE, masterRank,
                   tagData + 1, comm);
        } else if (rank == 2) {
          auto flattened = domainSolution.Flatten(top, 2);
          MPI_Send(flattened.data(), flattened.size(), MPI_DOUBLE, 3,
                   tagData + 2, comm);
        } else if (rank == 3) {
          auto size =
              ((domainSolution.GetM() + 1)) * ((domainSolution.GetN() + 1) - 2);
          Grid::line_t solution2(size);
          MPI_Recv(solution2.data(), solution2.size(), MPI_DOUBLE, 2,
                   tagData + 2, comm, status);
          auto halfPlaneSolution =
              domainSolution.Join(solution2, eDir::bottom, 2);
          auto flattened = halfPlaneSolution.Flatten(left, 2);
          MPI_Send(flattened.data(), flattened.size(), MPI_DOUBLE, masterRank,
                   tagData + 3, comm);
        }
        break;
      }
    }
    if (rank == masterRank) {
      auto size =
          ((domainSolution.GetM() + 1)) * ((domainSolution.GetN() + 1) - 2);
      Grid::line_t solVals(size);
      MPI_Recv(solVals.data(), solVals.size(), MPI_DOUBLE, nextRank,
               tagData + 1, comm, status);
      auto joinedGrid = domainSolution.Join(solVals, eDir::bottom, 2);

      size =
          (2 * (domainSolution.GetN() - 1) + 1) * (domainSolution.GetM() - 1);
      // size = 36;
      Grid::line_t rightHalf(size);
      MPI_Recv(rightHalf.data(), size, MPI_DOUBLE, prevRank, tagData + 3, comm,
               status);

      joinedGrid = joinedGrid.Join(rightHalf, eDir::right, 2);
      stop = MPI_Wtime();
      joinedGrid._execTime = std::chrono::duration_cast<Solution::time_t>(
          std::chrono::duration<double>(stop - start));
      joinedGrid.SaveToFile("mpi_4proc");
    }
  }
  MPI_Finalize();
  return 0;
}

/// @brief Calculate domain vertexes coordinates in case of 2 MPI process
std::array<int, 4> GetLimitsTwoProc(int rank, int M, int N) {
  // + 1 appears because of chosen splitting approach
  const int xMiddle = M / 2 + 1;
  int x0 = 0, xM = 0;
  switch (rank) {
  case 0:
    x0 = 0, xM = xMiddle;
    break;
  case 1:
    x0 = M - xMiddle, xM = M;
    break;
  }
  return {x0, xM, 0, N};
}

/// @brief Calculate domain vertexes coordinates in case of 4 MPI process
std::array<int, 4> GetLimitsFourProc(int rank, int M, int N) {
  // + 1 appears because of chosen splitting approach
  const int xMiddle = M / 2 + 1, yMiddle = N / 2 + 1;
  int x0 = 0, xM = 0, y0 = 0, yN = 0;
  switch (rank) {
  case 0:
    x0 = 0, xM = xMiddle, y0 = 0, yN = yMiddle;
    break;
  case 1:
    x0 = 0, xM = xMiddle, y0 = N - yMiddle, yN = N;
    break;
  case 2:
    x0 = M - xMiddle, xM = M, y0 = N - yMiddle, yN = N;
    break;
  case 3:
    x0 = M - xMiddle, xM = M, y0 = 0, yN = yMiddle;
    break;
  }
  return {x0, xM, y0, yN};
}

/// @brief Calculate domain coordinates, dividing original grid into pieces.
/// @param procNum Number of MPI processes (domains).
/// @param rank Current process id
/// @param M Original grid width
/// @param N Original grid height
/// @return Two vertexes: bottom left - (x0, y0), and top right - (xM, yN) - as
/// array [x0, xM, y0, yN]
std::array<int, 4> GetSectors(int procNum, int rank, int M, int N) {
  switch (procNum) {
  case 1:
    return {0, M, 0, N};
  case 2:
    return GetLimitsTwoProc(rank, M, N);
  case 4:
    return GetLimitsFourProc(rank, M, N);
  default:
    throw std::invalid_argument("Only 1, 2 or 4 processes could be launched\n");
  }
}