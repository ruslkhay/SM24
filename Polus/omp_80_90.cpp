#include "../src/inc.h"

int main() {

  //======OMP_NUM_THREADS======
  char *num_threads_str = getenv("OMP_NUM_THREADS");
  int num_threads = 2;

  if (num_threads_str != NULL) {
    num_threads = atoi(num_threads_str);
  }
  omp_set_dynamic(0);
  omp_set_num_threads(num_threads);
  std::cout << "Num threads: " << num_threads << std::endl;

  int maxIterations = 1e6;
  double tolerance = 1e-12;
  int x0 = 0, y0 = 0;
  int M = 80, N = 90;
  double h1 = 3.0 / M, h2 = 3.0 / N;
  auto s = Solution(M, N, x0, y0, h1, h2, maxIterations, tolerance);
  s.Find(omp, num_threads);
  s.SaveToFile("omp");
  return 0;
}