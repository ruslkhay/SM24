#include <filesystem>
#include <fstream>
#include <solution.hpp>

Grid::Grid(int M, int N) {
  _M = M;
  _N = N;
  _h1 =
      3.0 / M; // 3.0 is a number specific for my task (size of outer rectangle)
  _h2 = 3.0 / N;
  _nodes.assign(M, line_t(N, 0));
}

Grid::Grid(const matrix_t &grid) {
  _nodes = grid;
  _M = grid.size();
  _N = grid[0].size();
}

class Solution : public Grid {
  using time_t = std::chrono::microseconds;

private:
  std::filesystem::path _dirPath;
  int _maxIterations;
  time_t _execTime = time_t::duration::zero();
  double _tolerance;
  void CreateOutputDir(std::string buildDir = ".",
                       std::string outputDirName = "output");

public:
  Solution(int M, int N, int maxIterations, double tolerance);
  void SaveToFile(std::string fileName);
  void Find();
};

Solution::Solution(int M, int N, int maxIterations, double tolerance)
    : Grid(M + 1, N + 1) // +1 because if M = 3 => 4 nodes on X axis
{};

void Solution::Find() {
  auto start = std::chrono::high_resolution_clock::now();
  /* Calculations here */
  auto stop = std::chrono::high_resolution_clock::now();
  _execTime += std::chrono::duration_cast<time_t>(stop - start);
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
  fileName =
      fileName + "_" + std::to_string(_M) + '_' + std::to_string(_N) + ".txt";

  std::ofstream solutionFile(_dirPath / fileName);
  solutionFile << "Solution (" << _N << "," << _M << ")"
               << " Linear took " << _execTime.count() << " microsecs:\n\n";
  // Output the result
  for (int i = 0; i <= _M; i++) {
    for (int j = 0; j <= _N; j++) {
      solutionFile << std::fixed << std::setprecision(4) << _nodes[i][j] << " ";
    }
    solutionFile << std::endl;
  }
  solutionFile.close();
};
