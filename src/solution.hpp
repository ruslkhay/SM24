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

class Grid {

public:
  using line_t = std::vector<double>;
  using matrix_t = std::vector<line_t>;
  Grid(int M, int N);
  Grid(const matrix_t &grid);

  line_t GetColumn(int m);
  line_t GetRow(int n);

  line_t GetLeftBoarder() { return GetRow(0); };
  line_t GetRightBoarder() { return GetRow(_N - 1); };
  line_t GetTopBoarder() { return GetColumn(0); };
  line_t GetBottomBoarder() { return GetColumn(_M - 1); };

  inline matrix_t GetNodes() { return _nodes; }
  inline int GetM() { return _M; }
  inline int GetN() { return _N; }

  void SetLeftBoarder(const line_t &newBoarder);
  void SetRightBoarder(const line_t &newBoarder);
  void SetTopBoarder(const line_t &newBoarder);
  void SetBottomBoarder(const line_t &newBoarder);

  void Print();

protected:
  int _M;
  int _N;
  double _h1; // horizontal step
  double _h2; // vertical step
  matrix_t _nodes;
};

class Solution : public Grid {
public:
  using time_t = std::chrono::microseconds;

  Solution(int M, int N, int maxIterations, double tolerance);
  void SaveToFile(std::string fileName);
  void Find(sMethod method, int threads = 1);

private:
  int _maxIterations;
  time_t _execTime = time_t::duration::zero();
  double _tolerance;
  int _threads;

  matrix_t _a;
  matrix_t _b;
  matrix_t _F;

  std::filesystem::path _dirPath;
  void CreateOutputDir(std::string buildDir = ".",
                       std::string outputDirName = "output");
};

inline void Grid::SetTopBoarder(const line_t &newBoarder) {
  _nodes[0] = newBoarder;
}

inline void Grid::SetBottomBoarder(const line_t &newBoarder) {
  _nodes[_M] = newBoarder;
}

inline void Grid::SetLeftBoarder(const line_t &newBoarder) {
  for (int i = 0; i < _M + 1; ++i) {
    _nodes[i][0] = newBoarder[i];
  }
}

inline void Grid::SetRightBoarder(const line_t &newBoarder) {
  for (int i = 0; i < _M + 1; ++i) {
    _nodes[i][_N] = newBoarder[i];
  }
}

inline Grid::line_t Grid::GetColumn(int m) { return _nodes[m]; }

inline Grid::line_t Grid::GetRow(int n) {
  line_t res(_M);
  for (int i = 0; i < _M + 1; ++i) {
    res[i] = _nodes[i][n];
  }
  return res;
}

inline void Grid::Print() {
  std::string message;
  for (int i = 0; i < _M + 1; ++i) {
    for (int j = 0; j < _N + 1; ++j) {
      message += std::to_string(_nodes[i][j]) + "|";
    }
    message += '\n';
  }
  std::cout << message << std::endl;
}