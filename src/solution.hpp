#pragma once
#include "inc.h"

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
  inline int Size() { return (_M + 1) * (_N + 1); }

  void SetLeftBoarder(const line_t &newBoarder);
  void SetRightBoarder(const line_t &newBoarder);
  void SetTopBoarder(const line_t &newBoarder);
  void SetBottomBoarder(const line_t &newBoarder);

  void Print();
  line_t Flatten(eDir direction, int offset = 2);
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
  Solution Join(const line_t &boarderVal, eDir direction, int offset = 2);
  // Linear solution
  void ComputeABF() {
    ComputeA();
    ComputeB();
    ComputeF();
  };
  std::pair<double, double> CalculateTau();
  void CalculateResid();
  double CalculateMaxDiff(double tau);
  void OneStepOfSolution(double tau);

  // OpenMP solution
  void OMPComputeABF(int threads);
  std::pair<double, double> OMPCalculateTau();
  void OMPOneStepOfSolution(double tau);
  void OMPCalculateResid();
  double OMPCalculateMaxDiff(double tau);

  // Residuals manipulation

  matrix_t GetResiduals() { return _resid; }
  void SetResidColumn(int m, const line_t &valsNodesResid);
  void SetResidRow(int n, const line_t &valsNodesResid);
  line_t GetResidBoarder(eDir direction);
  void SetResidBoarder(eDir direction, const line_t &newBoarder);

  // Solution manipulation

  void SetSolutColumn(int m, const line_t &valsNodesResid);
  void SetSolutRow(int n, const line_t &valsNodesResid);
  line_t GetSolutBoarder(eDir direction);
  void SetSolutBoarder(eDir direction, const line_t &newBoarder);

  void Print();

  time_t _execTime = time_t::duration::zero();
  int _threads = 1;

private:
  int _maxIterations;
  double _tolerance;

  matrix_t _a;
  matrix_t _b;
  matrix_t _F;
  /// Residuals of approximation by chosen numeric schema
  matrix_t _resid;
  double _eps;

  void SetColumn(int m, matrix_t &mat, const line_t &newVals);
  void SetRow(int n, matrix_t &mat, const line_t &newVals);
  line_t GetColumn(int m, const matrix_t &vals);
  line_t GetRow(int n, const matrix_t &vals);

  void ComputeA();
  void ComputeB();
  void ComputeF();
  void ComputeW();
  void ComputeW(int procNum, int rank, int prevRank);

  void OMPComputeA();
  void OMPComputeB();
  void OMPComputeF();
  void OMPComputeW();

  // std::filesystem::path _dirPath;
  // void CreateOutputDir(std::string buildDir = ".",
  //                      std::string outputDirName = "output");
  double Product(const matrix_t &a, const matrix_t &b);
};

inline void Solution::SetColumn(int m, Grid::matrix_t &mat,
                                const Grid::line_t &newVals) {
  mat[m] = newVals;
};

inline Grid::line_t Solution::GetColumn(int m, const Grid::matrix_t &mat) {
  return mat[m];
};

inline void Solution::SetRow(int n, Grid::matrix_t &mat,
                             const Grid::line_t &newVals) {
  for (int i = 0; i < _M + 1; ++i) {
    mat[i][n] = newVals[i];
  }
};

inline Grid::line_t Solution::GetRow(int n, const Grid::matrix_t &mat) {
  line_t tmpBuf(_M + 1);
  for (int i = 0; i < _M + 1; ++i) {
    tmpBuf[i] = mat[i][n];
  }
  return tmpBuf;
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