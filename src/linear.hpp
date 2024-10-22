#pragma once
#include <vector>

namespace linear {
/// @brief Computer coefficient a_ij
/// @param a Matrix to store the values
void computeA(std::vector<std::vector<double>> &a);

void computeB(std::vector<std::vector<double>> &b);

void computeF(std::vector<std::vector<double>> &F);

inline double product(std::vector<std::vector<double>> &v1,
                      std::vector<std::vector<double>> &v2, double h1,
                      double h2);

void calculateW(const std::vector<std::vector<double>> &a,
                const std::vector<std::vector<double>> &b,
                const std::vector<std::vector<double>> &F,
                std::vector<std::vector<double>> &W, int maxIterations,
                double tolerance);
} // namespace linear
