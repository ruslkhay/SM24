#include <solution.hpp>

Grid::Grid(int M, int N, double h1, double h2) {
  _M = M;
  _N = N;
  _h1 = h1;
  _h2 = h2;
  _nodes.assign(M, std::vector<double>(N, 0));
}

Grid::Grid(int M, int N) {
  _M = M;
  _N = N;
  _nodes.assign(M, std::vector<double>(N, 0));
}

Grid::Grid(const std::vector<std::vector<double>> &grid) {
  _nodes = grid;
  _M = grid.size();
  _N = grid[0].size();
}

// class Solution
// {
// private:
//   int maxIterations = 1e5;
//   double tolerance = 1e-6;
// public:
//     Solution(/* args */);
//     ~Solution();
// };

// Solution::Solution(/* args */)
// {
// }

// Solution::~Solution()
// {
// }
