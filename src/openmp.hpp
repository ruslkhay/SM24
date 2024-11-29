#pragma once
#include "inc.h"

/// @brief Computer coefficient a_ij
/// @param a Matrix to store the values
/// @param TN Thread number to use
void computeA(std::vector<std::vector<double>> &a);

void computeB(std::vector<std::vector<double>> &b);

void computeF(std::vector<std::vector<double>> &F);

/// @brief Computer coefficients a_ij, b_ij, F_ij simultaneously.
/// It may seem odd, however it can boost performance while using threads:
/// call of unique computational function create threads each time, which is
/// greedy
/// @param a
/// @param b
/// @param F
/// @param TN Thread number
void computeJointABF(std::vector<std::vector<double>> &a,
                     std::vector<std::vector<double>> &b,
                     std::vector<std::vector<double>> &F);

double product(std::vector<std::vector<double>> &v1,
               std::vector<std::vector<double>> &v2, double h1, double h2);

void calculateW(const std::vector<std::vector<double>> &a,
                const std::vector<std::vector<double>> &b,
                const std::vector<std::vector<double>> &F,
                std::vector<std::vector<double>> &W, int maxIterations,
                double tolerance);
