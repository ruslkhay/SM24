#pragma once
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/// @brief Solution types. It's either linear (basic), OpenMP or MPI.
enum sMethod {
  lin,
  omp,
  mpi,
};

enum eDir {
  top,
  bottom,
  left,
  right,
};

class Grid {

public:
  using line_t = std::vector<double>;
  using matrix_t = std::vector<line_t>;
  Grid(int M, int N);
  Grid(int M, int N, double x0, double y0, double h1, double h2);
  Grid(const matrix_t &grid);

  line_t GetColumn(int m);
  line_t GetRow(int n);

  line_t GetTopBoarder() { return GetRow(0); };
  line_t GetBottomBoarder() { return GetRow(_N); };
  line_t GetLeftBoarder() { return GetColumn(0); };
  line_t GetRightBoarder() { return GetColumn(_M); };

  inline matrix_t GetNodes() { return _nodes; }
  inline int GetH1() { return _h1; }
  inline int GetH2() { return _h2; }
  inline int GetM() { return _M; }
  inline int GetN() { return _N; }

  void SetLeftBoarder(const line_t &newBoarder);
  void SetRightBoarder(const line_t &newBoarder);
  void SetTopBoarder(const line_t &newBoarder);
  void SetBottomBoarder(const line_t &newBoarder);

  void Print();
  /// @brief Convert 2D grid nodes into 1D vector
  /// @param direction Define what boarder values to get rid of.
  /// Dumping those are reasonable, because for Solution class all boarders are
  /// 0's
  line_t Flatten(eDir direction);
  // void Join(const line_t& flattened, eDir direction);
  // Grid Join(const line_t &flattened, eDir direction);
  Grid Join(const line_t &boarderVal, eDir direction);

protected:
  // number of columns
  int _M;
  /// number of rows
  int _N;
  int _x0, _y0;
  double _h1; // horizontal step
  double _h2; // vertical step
  matrix_t _nodes;
};

class Solution : public Grid {
public:
  using time_t = std::chrono::microseconds;

  Solution(int M, int N, double h1, double h2, double x0, double y0,
           int maxIterations, double tolerance);
  void SaveToFile(std::string fileName);
  void Find(sMethod method, int threads = 1);
  void ComputeA();
  void ComputeB();
  void ComputeF();
  Solution Join(const line_t &boarderVal, eDir direction);

private:
  int _maxIterations;
  time_t _execTime = time_t::duration::zero();
  double _tolerance;
  int _threads;

  matrix_t _a;
  matrix_t _b;
  matrix_t _F;
  double _eps;

  std::filesystem::path _dirPath;
  void CreateOutputDir(std::string buildDir = ".",
                       std::string outputDirName = "output");
};

inline void Grid::SetLeftBoarder(const line_t &newBoarder) {
  _nodes[0] = newBoarder;
}

inline void Grid::SetRightBoarder(const line_t &newBoarder) {
  _nodes[_M] = newBoarder;
}

inline void Grid::SetTopBoarder(const line_t &newBoarder) {
  for (int i = 0; i < _M + 1; ++i) {
    _nodes[i][0] = newBoarder[i];
  }
}

inline void Grid::SetBottomBoarder(const line_t &newBoarder) {
  for (int i = 0; i < _M + 1; ++i) {
    _nodes[i][_N] = newBoarder[i];
  }
}

inline Grid::line_t Grid::GetColumn(int m) { return _nodes[m]; }

inline Grid::line_t Grid::GetRow(int n) {
  line_t res(_M + 1);
  for (int i = 0; i < _M + 1; ++i) {
    res[i] = _nodes[i][n];
  }
  return res;
}

inline void Grid::Print() {
  std::string message;
  // Take specific index for each column
  for (int j = 0; j < _N + 1; ++j) {
    for (int i = 0; i < _M + 1; ++i) {
      message += std::to_string(_nodes[i][j]) + "|";
    }
    message += '\n';
  }
  std::cout << message << std::endl;
}
