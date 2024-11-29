#include "src/inc.h"

int main() {
  auto NMs = {std::pair<int, int>(10, 10), std::pair<int, int>(20, 20),
              std::pair<int, int>(40, 40)};
  int maxIterations = 1e5;
  int x0 = 0, y0 = 0;
  double tolerance = 1e-6;
  for (auto &[N, M] : NMs) {
    // Linear
    double h1 = 3.0 / M, h2 = 3.0 / N;
    auto s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
    s.Find(lin);
    s.SaveToFile("linear");
  }
  // OMP
  int M = 40, N = 40;
  double h1 = 3.0 / M, h2 = 3.0 / N;
  for (auto threadNum : {1, 4, 16}) {
    auto s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
    s.Find(omp, threadNum);
    s.SaveToFile("omp");
  };
  return 0;
}