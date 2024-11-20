#include <iostream>
#include <string>
#include <vector>

class Grid {
private:
  int _M;
  int _N;
  double _h1; // horizontal step
  double _h2; // vertical step
  std::vector<std::vector<double>> _nodes;

public:
  Grid(int M, int N);
  Grid(int M, int N, double h1, double h2);
  Grid(const std::vector<std::vector<double>> &grid);

  std::vector<double> GetLeftBoarder();
  std::vector<double> GetRightBoarder();
  std::vector<double> GetTopBoarder();
  std::vector<double> GetBottomBoarder();

  std::vector<double> GetColumn(int m);
  std::vector<double> GetRow(int n);

  std::vector<std::vector<double>> GetNodes();
  int GetM();
  int GetN();

  void SetLeftBoarder(const std::vector<double> &newBoarder);
  void SetRightBoarder(const std::vector<double> &newBoarder);
  void SetTopBoarder(const std::vector<double> &newBoarder);
  void SetBottomBoarder(const std::vector<double> &newBoarder);

  void Print();
};

inline void Grid::SetTopBoarder(const std::vector<double> &newBoarder) {
  _nodes[0] = newBoarder;
}

inline void Grid::SetBottomBoarder(const std::vector<double> &newBoarder) {
  _nodes[_M - 1] = newBoarder;
}

inline void Grid::SetLeftBoarder(const std::vector<double> &newBoarder) {
  for (int i = 0; i < _M; ++i) {
    _nodes[i][0] = newBoarder[i];
  }
}

inline void Grid::SetRightBoarder(const std::vector<double> &newBoarder) {
  for (int i = 0; i < _M; ++i) {
    _nodes[i][_N - 1] = newBoarder[i];
  }
}

inline std::vector<double> Grid::GetTopBoarder() { return _nodes[0]; }

inline std::vector<double> Grid::GetBottomBoarder() { return _nodes[_M - 1]; }

inline std::vector<double> Grid::GetLeftBoarder() {
  std::vector<double> res(_M);
  for (int i = 0; i < _M; ++i) {
    res[i] = _nodes[i][0];
  }
  return res;
}

inline std::vector<double> Grid::GetRightBoarder() {
  std::vector<double> res(_M);
  for (int i = 0; i < _M; ++i) {
    res[i] = _nodes[i][_N - 1];
  }
  return res;
}

inline std::vector<double> Grid::GetColumn(int m) { return _nodes[m]; }

inline std::vector<double> Grid::GetRow(int n) {
  std::vector<double> res(_M);
  for (int i = 0; i < _M; ++i) {
    res[i] = _nodes[i][n];
  }
  return res;
}

inline std::vector<std::vector<double>> Grid::GetNodes() { return _nodes; }
inline int Grid::GetM() { return _M; }
inline int Grid::GetN() { return _N; }

inline void Grid::Print() {
  std::string message;
  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      message += std::to_string(_nodes[i][j]) + "|";
    }
    message += '\n';
  }
  std::cout << message << std::endl;
}