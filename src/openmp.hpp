#pragma once
#include <vector>

/// @brief Computer coefficient a_ij
/// @param a Matrix to store the values
/// @param TN Thread number to use
void computeA(std::vector<std::vector<double>> &a, int TN);

void computeB(std::vector<std::vector<double>> &b, int TN);

void computeF(std::vector<std::vector<double>> &F, int TN);

inline double product(std::vector<std::vector<double>> &v1,
                      std::vector<std::vector<double>> &v2, double h1,
                      double h2, int TN = 1);

void calculateW(const std::vector<std::vector<double>> &a,
                const std::vector<std::vector<double>> &b,
                const std::vector<std::vector<double>> &F,
                std::vector<std::vector<double>> &W, int maxIterations,
                double tolerance, int TN);
