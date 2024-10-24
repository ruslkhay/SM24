#include "linear.hpp"
#include "auxility.hpp"
#include <cmath>
#include <iomanip> // for std::setprecision
#include <iostream>

void linear::computeA(std::vector<std::vector<double>> &a) {
  int M = a.size();
  int N = a[0].size() + 1;
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;
  double eps = std::pow(std::max(h1, h2), 2);
  // Calculate values for each inner node
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * h1, (j + 0.5) * h2};
      Point P_ij_next = {P_ij.x, P_ij.y + h2};

      double lv_ij = verticalShiftLen(P_ij, P_ij_next);
      // std::cout << "lv" << i+1 << j+1 << '=' << lv_ij << " h2=" << h2 <<
      // std::endl;
      if (lv_ij == h2) {
        a[i][j] = 1;
      } else {
        a[i][j] = lv_ij / h2 + (1 - lv_ij / h2) * (1 / eps);
      }
    }
  }
}

void linear::computeB(std::vector<std::vector<double>> &b) {
  int M = b.size() + 1;
  int N = b[0].size();
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;
  double eps = std::pow(std::max(h1, h2), 2);
  // Calculate values for each inner node
  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < N; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * h1, (j + 0.5) * h2};
      Point P_ji_next = {P_ij.x + h1, P_ij.y};

      double lh_ij = horizontalShiftLen(P_ij, P_ji_next);
      // std::cout << "lh" << i+1 << j+1 << '=' <<lh_ij << " h1=" << h1 << "
      // eps=" << eps << std::endl;
      if (lh_ij == h1) {
        b[i][j] = 1;
      } else {
        b[i][j] = lh_ij / h1 + (1 - lh_ij / h1) * (1 / eps);
      }
    }
  }
}

void linear::computeF(std::vector<std::vector<double>> &F) {
  int M = F.size() + 1;
  int N = F[0].size() + 1;
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;
  // Calculate values for each inner node
  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * h1, (j + 0.5) * h2};
      Point P_ij_diag = {P_ij.x + h1, P_ij.y + h2};

      F[i][j] = intersectionArea(P_ij, P_ij_diag) / (h1 * h2); // f(x*, y*) = 1
    }
  }
}

double linear::product(std::vector<std::vector<double>> &v1,
                       std::vector<std::vector<double>> &v2, double h1,
                       double h2) {
  double res = 0;
  for (int i = 0; i < static_cast<int>(v1.size()); i++) {
    for (int j = 0; j < static_cast<int>(v1[0].size()); j++) {
      res += h1 * h2 * v1[i][j] * v2[i][j];
    }
  }
  return res;
}

void linear::calculateW(const std::vector<std::vector<double>> &a,
                        const std::vector<std::vector<double>> &b,
                        const std::vector<std::vector<double>> &F,
                        std::vector<std::vector<double>> &W, int maxIterations,
                        double tolerance) {
  int M = F.size() + 1;
  int N = F[0].size() + 1;
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;

  std::vector<std::vector<double>> r(M + 1, std::vector<double>(N + 1, 0.0));
  std::vector<std::vector<double>> Ar(M + 1, std::vector<double>(N + 1, 0.0));
  // Perform the iterative steepest descent
  for (int iter = 0; iter < maxIterations; iter++) {
    std::vector<std::vector<double>> newW = W;
    double maxChange = 0.0;
    // Get residuals
    for (int i = 0; i < M - 1; i++) {
      for (int j = 0; j < N - 1; j++) {
        // Calculate the finite difference terms
        int I = i + 1;
        int J = j + 1;
        double term1 = (a[i][j] * (W[I][J] - W[I - 1][J])) / (h1 * h1);
        double term2 = (a[i + 1][j] * (W[I + 1][J] - W[I][J])) / (h1 * h1);
        double term3 = (b[i][j] * (W[I][J] - W[I][J - 1])) / (h2 * h2);
        double term4 = (b[i][j + 1] * (W[I][J + 1] - W[I][J])) / (h2 * h2);

        r[I][J] = (term2 - term1 + term4 - term3 + F[i][j]);
      }
    }
    // Get Ar
    for (int i = 0; i < M - 1; i++) {
      for (int j = 0; j < N - 1; j++) {
        int I = i + 1;
        int J = j + 1;
        // Calculate the finite difference terms
        double term1 = (a[i][j] * (r[I][J] - r[I - 1][J])) / (h1 * h1);
        double term2 = (a[i + 1][j] * (r[I + 1][J] - r[I][J])) / (h1 * h1);
        double term3 = (b[i][j] * (r[I][J] - r[I][J - 1])) / (h2 * h2);
        double term4 = (b[i][j + 1] * (r[I][J + 1] - r[I][J])) / (h2 * h2);

        Ar[I][J] = (term2 - term1 + term4 - term3 + F[i][j]);
      }
    }
    // Step of descend
    double theta = product(r, r, h1, h2) / product(Ar, r, h1, h2);
    for (int i = 0; i < M - 1; i++) {
      for (int j = 0; j < N - 1; j++) {
        int I = i + 1;
        int J = j + 1;
        newW[I][J] = W[I][J] - theta * r[I][J];
        // Track the maximum change
        maxChange = std::max(maxChange, std::abs(newW[i][j] - W[i][j]));
      }
    }
    // Update W after all computations
    W = newW;
    // Check for convergence
    if (maxChange < tolerance) {
      break;
    }
  }
}
