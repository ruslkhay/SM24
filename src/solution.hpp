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

  std::vector<double> GetColumn(int m);
  std::vector<double> GetRow(int n);

  std::vector<double> GetLeftBoarder() { return GetRow(0); };
  std::vector<double> GetRightBoarder() { return GetRow(_N - 1); };
  std::vector<double> GetTopBoarder() { return GetColumn(0); };
  std::vector<double> GetBottomBoarder() { return GetColumn(_M - 1); };

  inline std::vector<std::vector<double>> GetNodes() { return _nodes; }
  inline int GetM() { return _M; }
  inline int GetN() { return _N; }

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

inline std::vector<double> Grid::GetColumn(int m) { return _nodes[m]; }

inline std::vector<double> Grid::GetRow(int n) {
  std::vector<double> res(_M);
  for (int i = 0; i < _M; ++i) {
    res[i] = _nodes[i][n];
  }
  return res;
}

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