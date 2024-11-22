#include "solution.hpp"
#include "auxility.hpp"
#include "linear.hpp"
#include "omp.h"
#include "openmp.hpp"
#include <cmath>

Grid::Grid(int M, int N) {
  _M = M;
  _N = N;
  _h1 =
      3.0 / M; // 3.0 is a number specific for my task (size of outer rectangle)
  _h2 = 3.0 / N;
  _nodes.assign(_M + 1, line_t(_N + 1, 0));
}

Grid::Grid(const matrix_t &grid) {
  _nodes = grid;
  _M = grid.size() - 1;
  _N = grid[0].size() - 1;
}

// New
Grid::Grid(int M, int N, double x0, double y0, double h1, double h2) {
  _M = M, _N = N, _x0 = x0, _y0 = y0;
  _h1 = h1, _h2 = h2;
  _nodes.assign(_M + 1, line_t(_N + 1, 0));
}

Solution::Solution(int M, int N, double h1, double h2, double x0, double y0,
                   int maxIterations, double tolerance)
    : Grid(M, N, x0, y0, h1, h2) // +1 because if M = 3 => 4 nodes on X axis
{
  _maxIterations = maxIterations;
  _tolerance = tolerance;
  _eps = std::pow(std::max(_h1, _h2), 2);
  _a.assign(_M, line_t(_N - 1, 0));
  _b.assign(_M - 1, line_t(_N, 0));
  _F.assign(_M - 1, line_t(_N - 1, 0));
};

void Solution::Find(sMethod method, int threads) {
  _threads = threads;
  switch (method) {
  case lin: {
    auto start = std::chrono::high_resolution_clock::now();
    linear::computeA(_a);
    linear::computeB(_b);
    linear::computeF(_F);
    linear::calculateW(_a, _b, _F, _nodes, _maxIterations, _tolerance);
    auto stop = std::chrono::high_resolution_clock::now();
    _execTime = std::chrono::duration_cast<time_t>(stop - start);
  } break;
  case omp: {
    auto start = omp_get_wtime();
    computeJointABF(_a, _b, _F, threads);
    calculateW(_a, _b, _F, _nodes, _maxIterations, _tolerance, threads);
    auto stop = omp_get_wtime();
    _execTime = std::chrono::duration_cast<time_t>(
        std::chrono::duration<double>(stop - start));
  } break;
  case mpi:
    break;
  }
};

void Solution::CreateOutputDir(std::string buildDir,
                               std::string outputDirName) {
  const std::string buildFolder = buildDir;
  const std::string dirName = outputDirName;
  _dirPath = std::filesystem::path(buildFolder) / dirName;
  // Check if the directory already exists
  if (!std::filesystem::exists(_dirPath)) {
    // Create the directory
    if (std::filesystem::create_directory(_dirPath)) {
      std::cout << "Reports are stored in: " << _dirPath << std::endl;
    } else {
      std::cerr << "Failed to create directory: " << _dirPath << std::endl;
    }
  } else {
    std::cout << "Reports are stored in: " << _dirPath << std::endl;
  }
};

void Solution::SaveToFile(std::string fileName) {
  CreateOutputDir();
  fileName = fileName + "_" + std::to_string(_M) + '_' + std::to_string(_N) +
             '_' + std::to_string(_threads) + ".txt";

  std::ofstream solutionFile(_dirPath / fileName);
  solutionFile << "Solution (" << _N << "," << _M << ")"
               << " on " << _threads << " threads"
               << " took " << _execTime.count() << " microsecs:\n\n";
  // Output the result
  for (int i = 0; i <= _M; i++) {
    for (int j = 0; j <= _N; j++) {
      solutionFile << std::fixed << std::setprecision(4) << _nodes[i][j] << " ";
    }
    solutionFile << std::endl;
  }
  solutionFile.close();
};

void Solution::ComputeA() {
  // double h1 = 3.0 / M;
  // double h2 = 3.0 / N;
  // Calculate values for each inner node
  for (int i = 0; i < _M; i++) {
    for (int j = 0; j < _N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * _h1 + _x0, (j + 0.5) * _h2 + _y0};
      Point P_ij_next = {P_ij.x, P_ij.y + _h2};

      double lv_ij = verticalShiftLen(P_ij, P_ij_next);
      if (lv_ij == _h2) {
        _a[i][j] = 1;
      } else {
        _a[i][j] = lv_ij / _h2 + (1 - lv_ij / _h2) * (1 / _eps);
      }
    }
  }
}

void Solution::ComputeB() {
  // Calculate values for each inner node
  for (int i = 0; i < _M - 1; i++) {
    for (int j = 0; j < _N; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * _h1 + _x0, (j + 0.5) * _h2 + _y0};
      Point P_ji_next = {P_ij.x + _h1, P_ij.y};

      double lh_ij = horizontalShiftLen(P_ij, P_ji_next);
      if (lh_ij == _h1) {
        _b[i][j] = 1;
      } else {
        _b[i][j] = lh_ij / _h1 + (1 - lh_ij / _h1) * (1 / _eps);
      }
    }
  }
}

void Solution::ComputeF() {
  // Calculate values for each inner node
  for (int i = 0; i < _M - 1; i++) {
    for (int j = 0; j < _N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * _h1 + _x0, (j + 0.5) * _h2 + _y0};
      Point P_ij_diag = {P_ij.x + _h1, P_ij.y + _h2};

      _F[i][j] =
          intersectionArea(P_ij, P_ij_diag) / (_h1 * _h2); // f(x*, y*) = 1
    }
  }
}