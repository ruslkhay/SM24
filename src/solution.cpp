#include "solution.hpp"
#include "auxility.hpp"
#include "linear.hpp"
#include "omp.h"
#include "openmp.hpp"
#include <cmath>

// eDir GetOppositeDir(eDir direction){
//   eDir dir;
//   switch (direction)
//   {
//   case top:
//     dir = bottom;
//     break;
//   case bottom:
//     dir = top;
//     break;
//   case right:
//     dir = left;
//     break;
//   case left:
//     dir = right;
//     break;
//   }
//   return dir;
// };

Grid::Grid(int M, int N) {
  _M = M;
  _N = N;
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

/// @brief Convert 2D grid nodes into 1D vector
/// @param direction Define what boarder values to get rid of.
/// @param offset Define how many rows or columns needed to be dropped before
/// flattening Dumping those are reasonable, because for Solution class all
/// boarders are 0's
Grid::line_t Grid::Flatten(eDir direction, int offset) {
  std::vector<double> flattened;
  flattened.reserve((_M + 1) * (_N + 1));
  // Depending on the specified direction, remove the corresponding border
  switch (direction) {
  case left:
    for (int i = offset; i < _M + 1; i++) {
      for (int j = 0; j < _N + 1; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;

  case right:
    for (int i = 0; i < _M + 1 - offset; i++) {
      for (int j = 0; j < _N + 1; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;

  case top:
    for (int i = 0; i < _M + 1; i++) {
      for (int j = offset; j < _N + 1; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;

  case bottom:
    for (int i = 0; i < _M + 1; i++) {
      for (int j = 0; j < _N + 1 - offset; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;
  }

  return flattened;
}

Solution Solution::Join(const std::vector<double> &flattened, eDir direction) {
  int newM = _M; // Original columns, including borders
  int newN = _N; // Original rows, including borders
  int flatSize = flattened.size();

  if (flattened.size() % (_M + 1) != 0 && flattened.size() % (_N + 1) != 0) {
    throw std::invalid_argument(
        "The size of the flattened array does not match the expected size "
        "after border removal.");
  }
  int flatN, flatM;
  // Determine new dimensions based on the joining direction
  switch (direction) {
  case top:
  case bottom:
    --newN; // Get rid of top/bottom 0's boarder values of initial grid
    flatM = flatSize / _N; // TODO: add offset here <=> _N + 1 - offset
    flatN = flatSize / flatM;
    newN += flatN;
    break;
  case left:
  case right:
    --newM;
    flatN = flatSize / _M; // TODO: add offset here
    flatM = flatSize / flatN;
    newM += flatM;
    break;
  }
  // Grid joinedGrid(newM, newN);
  Solution joinedGrid(newM, newN, _x0, _y0, _h1, _h2, _maxIterations,
                      _tolerance);

  switch (direction) {
  case top: // Adding from the top
    for (int i = 0; i < newM + 1; ++i) {
      for (int j = 0; j < newN + 1; ++j) {
        if (j < flatN) {
          size_t flatIndex = i * flatN + j;
          joinedGrid._nodes[i][j] = flattened[flatIndex];
        } else {
          joinedGrid._nodes[i][j] = _nodes[i][j - flatN + 1];
        }
      }
    }
    break;

  case bottom: // Adding from the bottom
    for (int i = 0; i < newM + 1; ++i) {
      for (int j = 0; j < newN + 1; ++j) {
        if (j > newN - flatN) {
          size_t flatIndex = i * flatN + (j - (newN - flatN + 1));
          joinedGrid._nodes[i][j] = flattened[flatIndex];
        } else {
          joinedGrid._nodes[i][j] = _nodes[i][j];
        }
      }
    }
    break;

  case right: // Adding from the top
    for (int i = 0; i < newM + 1; ++i) {
      for (int j = 0; j < newN + 1; ++j) {
        if (i > newM - flatM) {
          size_t flatIndex = (i - (newM - flatM + 1)) * flatN + j;
          joinedGrid._nodes[i][j] = flattened[flatIndex];
        } else {
          joinedGrid._nodes[i][j] = _nodes[i][j];
        }
      }
    }
    break;

  case left: // Adding from the top
    for (int i = 0; i < newM + 1; ++i) {
      for (int j = 0; j < newN + 1; ++j) {
        if (i < flatM) {
          size_t flatIndex = i * flatN + j;
          joinedGrid._nodes[i][j] = flattened[flatIndex];
        } else {
          joinedGrid._nodes[i][j] = _nodes[i - flatM + 1][j];
        }
      }
    }
    break;
  }
  return joinedGrid;
}

Solution::Solution(int M, int N, double x0, double y0, double h1, double h2,
                   int maxIterations, double tolerance)
    : Grid(M, N, x0, y0, h1, h2) // +1 because if M = 3 => 4 nodes on X axis
{
  _maxIterations = maxIterations;
  _tolerance = tolerance;
  _eps = std::pow(std::max(_h1, _h2), 2);
  _a.assign(_M, line_t(_N - 1, 0));
  _b.assign(_M - 1, line_t(_N, 0));
  _F.assign(_M - 1, line_t(_N - 1, 0));
  _resid.assign(_M + 1, line_t(_N + 1, 0));
};

void Solution::Find(sMethod method, int threads) {
  _threads = threads;
  switch (method) {
  case lin: {
    auto start = std::chrono::high_resolution_clock::now();
    ComputeA();
    ComputeB();
    ComputeF();
    ComputeW();
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
  // Calculate values for each inner node
  for (int i = 0; i < _M; i++) {
    for (int j = 0; j < _N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5 + _x0) * _h1, (j + 0.5 + _y0) * _h2};
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
      Point P_ij = {(i + 0.5 + _x0) * _h1, (j + 0.5 + _y0) * _h2};
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
      Point P_ij = {(i + 0.5 + _x0) * _h1, (j + 0.5 + _y0) * _h2};
      Point P_ij_diag = {P_ij.x + _h1, P_ij.y + _h2};

      _F[i][j] =
          intersectionArea(P_ij, P_ij_diag) / (_h1 * _h2); // f(x*, y*) = 1
    }
  }
}

void Solution::CalculateResid() {
  for (int i = 0; i < _M - 1; i++) {
    for (int j = 0; j < _N - 1; j++) {
      int I = i + 1;
      int J = j + 1;
      double term1 =
          (_a[i][j] * (_nodes[I][J] - _nodes[I - 1][J])) / (_h1 * _h1);
      double term2 =
          (_a[i + 1][j] * (_nodes[I + 1][J] - _nodes[I][J])) / (_h1 * _h1);
      double term3 =
          (_b[i][j] * (_nodes[I][J] - _nodes[I][J - 1])) / (_h2 * _h2);
      double term4 =
          (_b[i][j + 1] * (_nodes[I][J + 1] - _nodes[I][J])) / (_h2 * _h2);

      _resid[I][J] = (term2 - term1 + term4 - term3 + _F[i][j]);
    }
  }
}

/// @brief Calculate step of descend for a numerical solution of the problem.
/// Initially calculate numerical schema of solution `Ar`. After that calculate
/// step `tau` of iterative descend
std::pair<double, double> Solution::CalculateTau() {
  matrix_t Ar(_M + 1, line_t(_N + 1, 0.0));
  CalculateResid();
  for (int i = 0; i < _M - 1; i++) {
    for (int j = 0; j < _N - 1; j++) {
      int I = i + 1;
      int J = j + 1;
      double term1 =
          (_a[i][j] * (_resid[I][J] - _resid[I - 1][J])) / (_h1 * _h1);
      double term2 =
          (_a[i + 1][j] * (_resid[I + 1][J] - _resid[I][J])) / (_h1 * _h1);
      double term3 =
          (_b[i][j] * (_resid[I][J] - _resid[I][J - 1])) / (_h2 * _h2);
      double term4 =
          (_b[i][j + 1] * (_resid[I][J + 1] - _resid[I][J])) / (_h2 * _h2);

      Ar[I][J] = (term2 - term1 + term4 - term3 + _F[i][j]);
    }
  }
  return {Product(_resid, _resid), Product(Ar, _resid)};
}

/// @brief Product in solution space
/// @param a first matrix
/// @param b second matrix
/// @attention parameters should be same size
double Solution::Product(const matrix_t &a, const matrix_t &b) {
  double res = 0;
  for (int i = 0; i < static_cast<int>(a.size()); i++) {
    for (int j = 0; j < static_cast<int>(a[0].size()); j++) {
      res += _h1 * _h2 * a[i][j] * b[i][j];
    }
  }
  return res;
}

/// @brief Make on step of iterative descent for problem solving
/// @return Norm of differences of solutions for neighboring steps
double Solution::OneStepOfSolution() {
  // Store differences between solution on different steps: w_(k+1) and w_k
  matrix_t diffs(_M + 1, line_t(_N + 1, 0.0));
  // Perform the iterative steepest descent
  auto [tau_nom, tau_denom] = CalculateTau();
  double tau = tau_nom / tau_denom;
  for (int i = 0; i < _M - 1; i++) {
    for (int j = 0; j < _N - 1; j++) {
      int I = i + 1;
      int J = j + 1;
      // Update values for solution
      _nodes[I][J] = _nodes[I][J] - tau * _resid[I][J];
      diffs[I][J] = tau * _resid[I][J];
    }
  }
  return std::sqrt(Product(diffs, diffs));
}

/// @brief Make on step of iterative descent for problem solving
/// @return Norm of differences of solutions for neighboring steps
double Solution::OneStepOfSolution(double tau) {
  // Store differences between solution on different steps: w_(k+1) and w_k
  matrix_t diffs(_M + 1, line_t(_N + 1, 0.0));
  // Perform the iterative steepest descent
  for (int i = 0; i < _M - 1; i++) {
    for (int j = 0; j < _N - 1; j++) {
      int I = i + 1;
      int J = j + 1;
      // Update values for solution
      _nodes[I][J] = _nodes[I][J] - tau * _resid[I][J];
      diffs[I][J] = tau * _resid[I][J];
    }
  }
  return std::sqrt(Product(diffs, diffs));
}

void Solution::ComputeW() {
  // Store differences between solution on different steps: w_(k+1) and w_k
  matrix_t diffs(_M + 1, line_t(_N + 1, 0.0));
  // Perform the iterative steepest descent
  for (int iter = 0; iter < _maxIterations; iter++) {
    double maxChange = 0.0;
    CalculateResid();
    auto [tau_nom, tau_denom] = CalculateTau();
    double tau = tau_nom / tau_denom;
    for (int i = 0; i < _M - 1; i++) {
      for (int j = 0; j < _N - 1; j++) {
        int I = i + 1;
        int J = j + 1;
        _nodes[I][J] = _nodes[I][J] - tau * _resid[I][J];
        diffs[I][J] = tau * _resid[I][J];
      }
    }
    maxChange = std::max(maxChange, std::sqrt(Product(diffs, diffs)));
    // Check for convergence
    if (maxChange < _tolerance) {
      break;
    }
  }
}