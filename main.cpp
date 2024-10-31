#include "src/auxility.hpp"
#include "src/linear.hpp"
#include "src/openmp.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>

int main() {

  // Specify the name of the directory to create in the build folder
  const std::string buildFolder =
      "."; // Replace with your build folder path if needed
  const std::string dirName =
      "output"; // Change this to your desired directory name
  const std::filesystem::path dirPath =
      std::filesystem::path(buildFolder) / dirName;

  // Check if the directory already exists
  if (!std::filesystem::exists(dirPath)) {
    // Create the directory
    if (std::filesystem::create_directory(dirPath)) {
      std::cout << "Reports are stored in: " << dirPath << std::endl;
    } else {
      std::cerr << "Failed to create directory: " << dirPath << std::endl;
    }
  } else {
    std::cout << "Reports are stored in: " << dirPath << std::endl;
  }
  auto NMs = {std::pair<int, int>(10, 10), std::pair<int, int>(20, 20),
              std::pair<int, int>(40, 40)};

  int maxIterations = 1e5;
  double tolerance = 1e-6;

  for (auto &[N, M] : NMs) {
    // Example matrices (a, b, F) initialized with some values
    std::vector<std::vector<double>> a(M, std::vector<double>(N - 1, 0.0));
    std::vector<std::vector<double>> b(M - 1, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> F(M - 1, std::vector<double>(N - 1, 0.0));
    std::vector<std::vector<double>> W(M + 1, std::vector<double>(N + 1, 0.0));

    std::string fileName =
        "linear_" + std::to_string(M) + '_' + std::to_string(N) + ".txt";

    std::ofstream solutionFile(dirPath / fileName);

    auto begin = std::chrono::high_resolution_clock::now();
    linear::computeA(a);
    linear::computeB(b);
    linear::computeF(F);
    linear::calculateW(a, b, F, W, maxIterations, tolerance);
    // After function call
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - begin);
    solutionFile << "Solution (" << N << "," << M << ")"
                 << " Linear took " << duration.count() << " microsecs:\n\n";
    // Output the result
    for (int i = 0; i <= M; i++) {
      for (int j = 0; j <= N; j++) {
        solutionFile << std::fixed << std::setprecision(4) << W[i][j] << " ";
      }
      solutionFile << std::endl;
    }
    solutionFile.close();

    // Back to default values
    for (int i = 0; i < M - 1; ++i) {
      std::fill(b[i].begin(), b[i].end(), 0.0);
      std::fill(F[i].begin(), F[i].end(), 0.0);
    }
    for (int i = 0; i < M; ++i) {
      std::fill(a[i].begin(), a[i].end(), 0.0);
    }
    for (int i = 0; i < M + 1; ++i) {
      std::fill(W[i].begin(), W[i].end(), 0.0);
    }
  }

  // // printABF(a, b, F);
  // // std::cout << "\n\n\n";

  int N = 40, M = 40;
  for (auto threadNum : {1, 4, 16}) {
    // Example matrices (a, b, F) initialized with some values
    std::vector<std::vector<double>> a(M, std::vector<double>(N - 1, 0.0));
    std::vector<std::vector<double>> b(M - 1, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> F(M - 1, std::vector<double>(N - 1, 0.0));
    std::vector<std::vector<double>> W(M + 1, std::vector<double>(N + 1, 0.0));

    auto fileName = "openmp_" + std::to_string(M) + '_' + std::to_string(N) +
                    "_TN" + std::to_string(threadNum) + ".txt";

    std::ofstream solutionFile(dirPath / fileName);

    auto start = omp_get_wtime();
    computeJointABF(a, b, F, threadNum);
    calculateW(a, b, F, W, maxIterations, tolerance, threadNum);
    auto end = omp_get_wtime();
    // Output the result
    solutionFile << std::fixed << "Solution (" << N << "," << M << ")"
                 << " OpenMP took (" << threadNum << " threads) "
                 << static_cast<uint>((end - start) * 1e6) << " microsecs:\n\n";
    for (int i = 0; i <= M; i++) {
      for (int j = 0; j <= N; j++) {
        solutionFile << std::fixed << std::setprecision(4) << W[i][j] << " ";
      }
      solutionFile << std::endl;
    }
    solutionFile.close();

    // Back to default values
    for (int i = 0; i < M - 1; ++i) {
      std::fill(b[i].begin(), b[i].end(), 0.0);
      std::fill(F[i].begin(), F[i].end(), 0.0);
    }
    for (int i = 0; i < M; ++i) {
      std::fill(a[i].begin(), a[i].end(), 0.0);
    }
    for (int i = 0; i < M + 1; ++i) {
      std::fill(W[i].begin(), W[i].end(), 0.0);
    }
  }
  return 0;
}