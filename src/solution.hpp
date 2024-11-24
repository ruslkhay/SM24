#pragma once
#include <chrono>
#include <cmath>
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

  line_t GetTopBoarder() { return GetRow(1); };
  line_t GetBottomBoarder() { return GetRow(_N - 1); };
  line_t GetLeftBoarder() { return GetColumn(1); };
  line_t GetRightBoarder() { return GetColumn(_M - 1); };

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
  line_t Flatten(eDir direction, int offset = 1);
  Grid Join(const line_t &boarderVal, eDir direction);

protected:
  int _M;       /// Number of columns
  int _N;       /// Number of rows
  int _x0, _y0; /// Bottom left corner position for grid
  double _h1;   /// Horizontal step of grid
  double _h2;   /// Vertical step of grid
  /// Values, stored in grid nodes.
  /// In case of `Solution` class, it is solution itself
  matrix_t _nodes;
};

class Solution : public Grid {
public:
  using time_t = std::chrono::microseconds;

  Solution(int M, int N, double h1, double h2, double x0, double y0,
           int maxIterations, double tolerance);
  void SaveToFile(std::string fileName);
  void Find(sMethod method, int threads = 1);

  line_t GetTopResid();
  line_t GetBottomResid();
  line_t GetLeftResid();
  line_t GetRightResid();

  void SetLeftResid(const line_t &newBoarder);
  void SetRightResid(const line_t &newBoarder);
  void SetTopResid(const line_t &newBoarder);
  void SetBottomResid(const line_t &newBoarder);

  Solution Join(const line_t &boarderVal, eDir direction, int offset = 1);
  void ComputeABF() {
    ComputeA();
    ComputeB();
    ComputeF();
  };
  void CalculateResid();
  std::pair<double, double> CalculateTau();
  matrix_t GetResiduals() { return _resid; }
  double OneStepOfSolution(double tau);
  double OneStepOfSolution();

private:
  int _maxIterations;
  time_t _execTime = time_t::duration::zero();
  double _tolerance;
  int _threads = 1;

  matrix_t _a;
  matrix_t _b;
  matrix_t _F;
  /// Residuals of approximation by chosen numeric schema
  matrix_t _resid;
  double _eps;

  void ComputeA();
  void ComputeB();
  void ComputeF();
  void ComputeW();
  void ComputeW(int procNum, int rank, int prevRank);

  std::filesystem::path _dirPath;
  void CreateOutputDir(std::string buildDir = ".",
                       std::string outputDirName = "output");
  double Product(const matrix_t &a, const matrix_t &b);
};

inline Solution::line_t Solution::GetLeftResid() { return _resid[1]; }

inline Solution::line_t Solution::GetRightResid() { return _resid[_M - 1]; }

inline Solution::line_t Solution::GetTopResid() {
  line_t res(_M + 1);
  for (int i = 0; i < _M + 1; ++i) {
    res[i] = _resid[i][1];
  }
  return res;
}

inline Solution::line_t Solution::GetBottomResid() {
  line_t res(_M + 1);
  for (int i = 0; i < _M + 1; ++i) {
    res[i] = _resid[i][_N - 1];
  }
  return res;
}

inline void Solution::SetLeftResid(const line_t &newBoarder) {
  _resid[0] = newBoarder;
}

inline void Solution::SetRightResid(const line_t &newBoarder) {
  _resid[_M] = newBoarder;
}

inline void Solution::SetTopResid(const line_t &newBoarder) {
  for (int i = 0; i < _M + 1; ++i) {
    _resid[i][0] = newBoarder[i];
  }
}

inline void Solution::SetBottomResid(const line_t &newBoarder) {
  for (int i = 0; i < _M + 1; ++i) {
    _resid[i][_N] = newBoarder[i];
  }
}

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
