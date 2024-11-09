#pragma once
#include <vector>

namespace linear {
/// @brief Computer coefficient a_ij
/// @param a Matrix to store the values
void computeA(std::vector<std::vector<double>> &a);

/// @brief Computer coefficient b_ij
/// @param b Matrix to store the values
void computeB(std::vector<std::vector<double>> &b);

/// @brief Computer coefficient F_ij
/// @param F Matrix to store the values
void computeF(std::vector<std::vector<double>> &F);

/// @brief Dot product in vector space described in difference schema
/// @param v1 First vector
/// @param v2 Second vector
/// @param h1 Grid step in horizontal direction
/// @param h2 Grid step in vertical direction
/// @return Dot product of two vectors defined on gird with known parameters
inline double product(std::vector<std::vector<double>> &v1,
                      std::vector<std::vector<double>> &v2, double h1,
                      double h2);

/// @brief Computer solution w_ij
/// @param a Matrix to store the values a_ij
/// @param b Matrix to store the values b_ij
/// @param F Matrix to store the values F_ij
/// @param W Matrix to store the solution values w_ij
/// @param maxIterations Limit on iterations while building solution
/// @param tolerance Limit on accuracy of solution
void calculateW(const std::vector<std::vector<double>> &a,
                const std::vector<std::vector<double>> &b,
                const std::vector<std::vector<double>> &F,
                std::vector<std::vector<double>> &W, int maxIterations,
                double tolerance);
} // namespace linear
