#include "solution.hpp"
#include "auxility.hpp"
#include "linear.hpp"
#include "omp.h"
#include "openmp.hpp"
#include <cmath>

Grid::Grid(int M, int N) {
  _M = M;
  _N = N;
  // _h1 =
  // 3.0 / M; // 3.0 is a number specific for my task (size of outer rectangle)
  // _h2 = 3.0 / N;
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

Grid::line_t Grid::Flatten(eDir direction) {
  // Create a vector to hold the flattened results
  std::vector<double> flattened;
  flattened.reserve(_M * _N);

  // Depending on the specified direction, remove the corresponding border
  switch (direction) {
  case left: // Remove left border
    for (int i = 1; i < _M + 1; i++) {
      for (int j = 0; j < _N + 1; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;

  case right: // Remove right border
    for (int i = 0; i < _M; i++) {
      for (int j = 0; j < _N + 1; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;

  case top: // Remove top border
    for (int i = 0; i < _M + 1; i++) {
      for (int j = 1; j < _N + 1; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;

  case bottom: // Remove bottom border
    for (int i = 0; i < _M + 1; i++) {
      for (int j = 0; j < _N; j++) {
        flattened.push_back(_nodes[i][j]);
      }
    }
    break;
  }

  return flattened;
}

// bool shouldFillOriginal(int i, int j, eDir direction, int flatSize, int _N,
// int _M) {
//     // Logic to determine if we should fill original grid values based on the
//     direction switch (direction) {
//         case top: return i >= (flatSize / _N);
//         case bottom: return true; // Always fill original values if going to
//         bottom case right: return j < _N; case left: return j >= (flatSize /
//         _M);
//     }
//     return false;
// }

// bool isInsideOriginalGrid(int i, int j, eDir direction, int flatSize, int _N,
// int _M) {
//     // Logic to check if the indices are inside the original grid after
//     offset switch (direction) {
//         case top: return i >= (flatSize / _N);
//         case bottom: return i < _M;
//         case right: return j < _N;
//         case left: return j < _N;
//     }
//     return false;
// }

// int offsetRow(int i, eDir direction, int flatSize, int _N, int _M) {
//     // Calculate row offset depending on joining direction
//     switch (direction) {
//         case top: return flatSize / _N;
//         case bottom: return 0;
//         case right: return 0;
//         case left: return 0;
//     }
//     return 0;
// }

// int offsetCol(int j, eDir direction, int flatSize, int _N, int _M) {
//     // Calculate column offset depending on joining direction
//     switch (direction) {
//         case top: return 0;
//         case bottom: return 0;
//         case right: return 0;
//         case left: return flatSize / _M;
//     }
//     return 0;
// }

// Grid Grid::Join(const std::vector<double>& flattened, eDir direction) {
//   int newRows = _M;
//   int newCols = _N;
//   int flatSize = flattened.size();

//   // Determine new dimensions based on the joining direction
//   switch (direction) {
//     case top: // Top
//       newRows += (flatSize / _N);
//       break;
//     case bottom: // Bottom
//       newRows += (flatSize / _N);
//       break;
//     case right: // Right
//       newCols += (flatSize / _M);
//       break;
//     case left: // Left
//       newCols += (flatSize / _M);
//       break;
//   }

//   // Create a new grid with the updated dimensions
//   Grid newGrid(newRows, newCols);

//   // Fill the new grid with appropriate values
//   int flatIndex = 0;

//   printf("newN = %d, newM = %d\n", newCols, newRows);
//   for (int i = 0; i < newRows; ++i) {
//     for (int j = 0; j < newCols; ++j) {
//       if (shouldFillOriginal(i, j, direction, flatSize, _N,_M)) {
//         // Fill from the original grid
//         if (isInsideOriginalGrid(i, j, direction, flatSize, _N,_M)) {
//           int newI = i - offsetRow(i, direction, flatSize, _N,_M);
//           int newJ = j - offsetCol(j, direction, flatSize, _N,_M);
//           newGrid._nodes[i][j] = _nodes[newI][newJ];
//         } else {
//           newGrid._nodes[i][j] = 7; // Fill with 7 if outside original grid
//         }
//       } else {
//         // Fill with flattened values if there's space
//         if (flatIndex < flatSize) {
//           newGrid._nodes[i][j] = flattened[flatIndex++];
//         } else {
//           newGrid._nodes[i][j] = 7; // Fill with 7 if not enough flattened
//           values
//         }
//       }
//     }
//   }
//   return newGrid;
// }

Grid Grid::Join(const std::vector<double> &flattened, eDir direction) {
  int newM = _M; // Original columns, including borders
  int newN = _N; // Original rows, including borders
  int flatSize = flattened.size();

  if (flattened.size() % _M != 0 && flattened.size() % _N != 0) {
    throw std::invalid_argument(
        "The size of the flattened array does not match the expected size "
        "after border removal.");
  }
  int flatN = 0, flatM = 0;
  // Determine new dimensions based on the joining direction
  switch (direction) {
  case top:
  case bottom:
    --newN; // Get rid of top/bottom 0's boarder values of initial grid
    flatM = flatSize / _N;
    flatN = flatSize / flatM;
    newN += flatN;
    break;
  case left:
  case right:
    --newM;
    flatN = flatSize / _M;
    flatM = flatSize / flatN;
    newM += flatM;
    break;
  }
  Grid joinedGrid(newM, newN);

  printf("newM=%d, newN=%d, flatSize=%d, flatM=%d, flatN=%d\n", newM, newN,
         flatSize, flatM, flatN);
  switch (direction) {
  case top: // Adding from the top
    for (int i = 0; i < newM + 1; ++i) {
      for (int j = 0; j < newN + 1; ++j) {
        // printf("i=%d, j=%d\n", i, j);
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
          // printf("flatIndex = %ld, flatten=%f\n", flatIndex,
          //  flattened[flatIndex]);
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
          // printf("flatIndex = %ld, flatten=%f\n", flatIndex,
          //  flattened[flatIndex]);
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
          // printf("flatIndex = %ld, flatten=%f\n", flatIndex,
          //  flattened[flatIndex]);
          joinedGrid._nodes[i][j] = 10 * flattened[flatIndex];
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
};

void Solution::Find(sMethod method, int threads) {
  _threads = threads;
  switch (method) {
  case lin: {
    auto start = std::chrono::high_resolution_clock::now();
    // linear::computeA(_a);
    ComputeA();
    ComputeB();
    // linear::computeB(_b);
    ComputeF();
    // linear::computeF(_F);
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