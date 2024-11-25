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

class Solution : public Grid {
public:
  using time_t = std::chrono::microseconds;

  Solution(int M, int N, double h1, double h2, double x0, double y0,
           int maxIterations, double tolerance);
  void SaveToFile(std::string fileName);
  void Find(sMethod method, int threads = 1);
  Solution Join(const line_t &boarderVal, eDir direction, int offset = 1);
  void ComputeABF() {
    ComputeA();
    ComputeB();
    ComputeF();
  };
  matrix_t GetResiduals() { return _resid; }
  std::pair<double, double> CalculateTau();
  double OneStepOfSolution(double tau);
  double OneStepOfSolution();

  std::pair<line_t, line_t> GetColumn(int m);
  std::pair<line_t, line_t> GetRow(int n);

  void SetColumn(int m, const std::pair<line_t, line_t> &valsNodesResid);
  void SetRow(int n, const std::pair<line_t, line_t> &valsNodesResid);

  // We are interested only in getting inner nodes
  std::pair<line_t, line_t> GetTopBoarder() { return GetRow(1); };
  std::pair<line_t, line_t> GetBottomBoarder() { return GetRow(_N - 1); };
  std::pair<line_t, line_t> GetLeftBoarder() { return GetColumn(1); };
  std::pair<line_t, line_t> GetRightBoarder() { return GetColumn(_M - 1); };
  // We are interested only in changing outer nodes boarder values
  void SetTopBoarder(const std::pair<line_t, line_t> &newBoarder) {
    SetRow(0, newBoarder);
  };
  void SetBottomBoarder(const std::pair<line_t, line_t> &newBoarder) {
    SetRow(_N, newBoarder);
  };
  void SetLeftBoarder(const std::pair<line_t, line_t> &newBoarder) {
    SetColumn(0, newBoarder);
  };
  void SetRightBoarder(const std::pair<line_t, line_t> &newBoarder) {
    SetColumn(_M, newBoarder);
  };

  void Print();

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
  void CalculateResid();
  double Product(const matrix_t &a, const matrix_t &b);
};

inline std::pair<Grid::line_t, Grid::line_t> Solution::GetColumn(int m) {
  return {_nodes[m], _resid[m]};
};

inline std::pair<Grid::line_t, Grid::line_t> Solution::GetRow(int n) {
  line_t valsNodes(_M + 1);
  line_t valsResid(_M + 1);
  for (int i = 0; i < _M + 1; ++i) {
    valsNodes[i] = _nodes[i][n];
    valsResid[i] = _nodes[i][n];
  }
  return {valsNodes, valsResid};
};

inline void Solution::SetColumn(
    int m, const std::pair<Grid::line_t, Grid::line_t> &valsNodesResid) {
  _nodes[m] = valsNodesResid.first;
  _resid[m] = valsNodesResid.second;
};

inline void
Solution::SetRow(int n,
                 const std::pair<Grid::line_t, Grid::line_t> &valsNodesResid) {
  for (int i = 0; i < _M + 1; ++i) {
    _nodes[i][n] = valsNodesResid.first[i];
    _resid[i][n] = valsNodesResid.second[i];
  }
};

inline void Solution::Print() {
  std::string message = "\tSolution:\n";
  // take specific index for each column
  for (int j = 0; j < _N + 1; ++j) {
    for (int i = 0; i < _M + 1; ++i) {
      message += "\t" + std::to_string(_nodes[i][j]) + "|";
    }
    message += '\n';
  }
  message += "\tResiduals:\n";
  for (int j = 0; j < _N + 1; ++j) {
    for (int i = 0; i < _M + 1; ++i) {
      message += "\t" + std::to_string(_resid[i][j]) + "|";
    }
    message += '\n';
  }
  std::cout << message << std::endl;
}