#include "src/linear.hpp"
#include "src/openmp.hpp"
#include <benchmark/benchmark.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>

static void BM_CalculateW_OpenMP_40_40(benchmark::State &state) {
  int M = 40;
  int N = 40;
  int threadNum = state.range(0);

  std::vector<std::vector<double>> a(M, std::vector<double>(N - 1, 0.0));
  std::vector<std::vector<double>> b(M - 1, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> F(M - 1, std::vector<double>(N - 1, 0.0));
  std::vector<std::vector<double>> W(M + 1, std::vector<double>(N + 1, 0.0));
  int maxIterations = 1000;
  double tolerance = 1e-5;

  for (auto _ : state) {
    computeA(a, threadNum);
    computeB(b, threadNum);
    computeF(F, threadNum);
    calculateW(a, b, F, W, maxIterations, tolerance, threadNum);
    // state.SetIterationTime(time);
  }
}
BENCHMARK(BM_CalculateW_OpenMP_40_40)->Arg(1)->Arg(4)->Arg(16)->UseRealTime();

static void BM_CalculateW_Linear(benchmark::State &state) {
  int M = state.range(0);
  int N = state.range(1);

  std::vector<std::vector<double>> a(M, std::vector<double>(N - 1, 0.0));
  std::vector<std::vector<double>> b(M - 1, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> F(M - 1, std::vector<double>(N - 1, 0.0));
  std::vector<std::vector<double>> W(M + 1, std::vector<double>(N + 1, 0.0));
  int maxIterations = 1000;
  double tolerance = 1e-5;

  for (auto _ : state) {
    linear::computeA(a);
    linear::computeB(b);
    linear::computeF(F);
    linear::calculateW(a, b, F, W, maxIterations, tolerance);
  }
}
BENCHMARK(BM_CalculateW_Linear)
    ->Args({10, 10})
    ->Args({20, 20})
    ->Args({40, 40})
    ->UseRealTime();

// Register the benchmark
BENCHMARK_MAIN();