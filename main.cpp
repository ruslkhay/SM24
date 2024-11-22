#include "src/auxility.hpp"
#include "src/linear.hpp"
#include "src/openmp.hpp"
#include "src/solution.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>

int main() {

  auto NMs = {std::pair<int, int>(10, 10), std::pair<int, int>(20, 20),
              std::pair<int, int>(40, 40)};

  int maxIterations = 1e5;
  double tolerance = 1e-6;
  int x0 = 0, y0 = 0;
  auto meth = sMethod::lin;
  for (auto &[N, M] : NMs) {
    double h1 = 3.0 / M, h2 = 3.0 / N;
    auto s = Solution(M, N, h1, h2, x0, y0, maxIterations, tolerance);
    s.Find(meth);
    s.SaveToFile("linear");
  }

  int M = 40, N = 40;
  double h1 = 3.0 / M, h2 = 3.0 / N;
  meth = sMethod::omp;
  for (auto threadNum : {1, 4, 16}) {
    auto s = Solution(M, N, h1, h2, x0, y0, maxIterations, tolerance);
    // auto s = Solution(M, N, maxIterations, tolerance);
    s.Find(meth, threadNum);
    s.SaveToFile("openmp");
  }
  return 0;
}