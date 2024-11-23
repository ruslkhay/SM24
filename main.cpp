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

  auto NMs = {std::pair<int, int>(4, 4), std::pair<int, int>(20, 20),
              std::pair<int, int>(40, 40)};

  int maxIterations = 1e5;
  double tolerance = 1e-6;
  int x0 = 0, y0 = 0;
  auto meth = sMethod::lin;
  for (auto &[N, M] : NMs) {
    double h1 = 3.0 / M, h2 = 3.0 / N;
    auto s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
    s.Find(meth);
    // s.SaveToFile("linear");
    std::cout << "Initial grid:\n";
    s.Print();
    std::cout << "Flattened grid:\n";
    auto flattened = s.Flatten(eDir::right);
    for (auto elem : flattened) {
      std::cout << elem << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    auto r = s.Join(flattened, eDir::left);
    std::cout << "Joined grid:\n";
    r.Print();
    break;
  }

  // int M = 40, N = 40;
  // double h1 = 3.0 / M, h2 = 3.0 / N;
  // meth = sMethod::omp;
  // for (auto threadNum : {1, 4, 16}) {
  //   auto s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
  //   // auto s = Solution(M, N, maxIterations, tolerance);
  //   s.Find(meth, threadNum);
  //   s.SaveToFile("openmp");
  // }
  return 0;
}