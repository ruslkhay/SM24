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

  int maxIterations = 1e5;
  int x0 = 0, y0 = 0;
  // for (auto &[N, M] : NMs) {
  for (int i = 1; i <= 3; ++i) {
    // Linear
    int M = 40 * i, N = 40 * i;
    // double tolerance =  1e-6;
    double tolerance = std::pow(10, -(5 + i));
    double h1 = 3.0 / M, h2 = 3.0 / N;
    auto s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
    s.Find(lin);
    s.SaveToFile("linear");
    // // 1 OMP thread (as linear)
    // s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
    // s.Find(omp, 1);
    // s.SaveToFile("omp");
    // // OMP
    // for (auto threadNum : {2, 4, 8, 16}) {
    //   s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
    //   s.Find(omp, threadNum * i);
    //   s.SaveToFile("omp");
    // };
  }
  return 0;
}