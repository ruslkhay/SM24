#include <cmath>
#include <iomanip> // for std::setprecision
#include <iostream>
#include <omp.h> // Include the OpenMP header
#include <vector>

struct Point {
  double x, y;
};

bool isInTrapezoid(Point p) {
  return (p.y <= -3 * p.x + 9 && p.x >= 2 && p.x <= 3) ||
         (p.y >= 0 && p.y <= 3 && p.x >= 0 && p.x <= 2);
}

bool isInRectangle(Point p) {
  return p.y >= 0 && p.y <= 3 && p.x >= 0 && p.x <= 3;
}

// BC side lays on y = -3x + 9 line
double lineBC(char axis, double val) {
  if (axis == 'x') {
    return -3 * val + 9;
  } else {
    return 3 - val / 3;
  }
}

double horizontalShiftLen(Point left, Point right) {
  if (isInTrapezoid(right)) {
    return right.x - left.x;
  } else {
    auto rLimit = lineBC('y', right.y);
    return rLimit - left.x;
  }
}

double verticalShiftLen(Point bottom, Point top) {
  if (isInTrapezoid(top)) {
    return top.y - bottom.y;
  } else {
    auto tLimit = lineBC('x', top.x);
    return tLimit - bottom.y;
  }
}

void computeA(std::vector<std::vector<double>> &a) {
  int M = a.size();
  int N = a[0].size() + 1;
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;
  double eps = std::pow(std::max(h1, h2), 2);

// Parallelizing the outer loop with OpenMP
#pragma omp parallel for num_threads(1)
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * h1, (j + 0.5) * h2};
      Point P_ij_next = {P_ij.x, P_ij.y + h2};

      double lv_ij = verticalShiftLen(P_ij, P_ij_next);
      if (lv_ij == h2) {
        a[i][j] = 1;
      } else {
        a[i][j] = lv_ij / h2 + (1 - lv_ij / h2) * (1 / eps);
      }
    }
  }
}

void computeB(std::vector<std::vector<double>> &b) {
  int M = b.size() + 1;
  int N = b[0].size();
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;
  double eps = std::pow(std::max(h1, h2), 2);

// Parallelizing the outer loop with OpenMP
#pragma omp parallel for num_threads(1)
  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < N; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * h1, (j + 0.5) * h2};
      Point P_ji_next = {P_ij.x + h1, P_ij.y};

      double lh_ij = horizontalShiftLen(P_ij, P_ji_next);
      if (lh_ij == h1) {
        b[i][j] = 1;
      } else {
        b[i][j] = lh_ij / h1 + (1 - lh_ij / h1) * (1 / eps);
      }
    }
  }
}

/// @brief Calculate area for the gird rectangles, which is set by adjacent
double intersectionArea(Point bottomLeftPoint, Point topRightPoint) {
  double h1 = topRightPoint.x - bottomLeftPoint.x;
  double h2 = topRightPoint.y - bottomLeftPoint.y;
  double S_ij = 0; // f(x*, y*) = 1
  if (isInTrapezoid(topRightPoint)) {
    S_ij = h1 * h2;
  } else if (!isInTrapezoid(topRightPoint) && isInTrapezoid(bottomLeftPoint) &&
             (bottomLeftPoint.y != -3 * bottomLeftPoint.x + 9)) {
    // If only one point is inside the trapezoid, calculate the distance to the
    // edge
    Point topLeftPoint = {bottomLeftPoint.x, bottomLeftPoint.y + h2};
    Point bottomRightPoint = {bottomLeftPoint.x + h1, bottomLeftPoint.y};

    Point p1, p2; // p1.y > p2.y, p1.x < p2.x always
    if (isInTrapezoid(topLeftPoint)) {
      p1 = {lineBC('y', topRightPoint.y), topLeftPoint.y};
    } else {
      p1 = {topLeftPoint.x, lineBC('x', topLeftPoint.x)};
    }
    if (isInTrapezoid(bottomRightPoint)) {
      p2 = {bottomRightPoint.x, lineBC('x', bottomRightPoint.x)};
    } else {
      p2 = {lineBC('y', bottomRightPoint.y), bottomRightPoint.y};
    }
    bool p1OnVert = (p1.y < topLeftPoint.y && p1.x == topLeftPoint.x);
    bool p1OnHor = (p1.y == topLeftPoint.y && p1.x > topLeftPoint.x);
    bool p1InCorner = (p1.x == topLeftPoint.x && p1.y == topLeftPoint.y);
    bool p2OnVert = (p2.y > bottomRightPoint.y && p2.x == bottomRightPoint.x);
    bool p2OnHor = (p2.y == bottomRightPoint.y && p2.x < bottomRightPoint.x);
    bool p2InCorner =
        (p2.x == bottomRightPoint.x && p2.y == bottomRightPoint.y);

    // If interception is a triangle
    if ((p1OnVert | p1InCorner) && (p2OnHor | p2InCorner)) {
      auto sideV = p1.y - bottomLeftPoint.y;
      auto sideH = p2.x - bottomLeftPoint.x;
      S_ij = sideV * sideH * 1 / 2;
    }
    // If grid rectangle without interception (complement of interception) is a
    // triangle
    if ((p1OnHor | p1InCorner) && (p2OnVert | p2InCorner)) {
      auto sideV = topRightPoint.x - p1.x;
      auto sideH = topRightPoint.y - p2.y;
      S_ij = h1 * h2 - sideV * sideH * 1 / 2;
    }
    // If S_ij is a trapezoid with vertical bases
    if (p1OnVert && p2OnVert) {
      auto h = bottomRightPoint.x - bottomLeftPoint.x;
      auto base1 = p2.y - bottomRightPoint.y;
      auto base2 = p1.y - bottomLeftPoint.y;
      S_ij = h * (base1 + base2) / 2;
    }
    // If S_ij is a trapezoid with horizontal bases
    if (p1OnHor && p2OnHor) {
      auto h = topLeftPoint.y - bottomLeftPoint.y;
      auto base1 = p1.x - topLeftPoint.x;
      auto base2 = p2.x - bottomLeftPoint.x;
      S_ij = h * (base1 + base2) / 2;
    }
  }
  return S_ij / (h1 * h2);
}

void computeF(std::vector<std::vector<double>> &F) {
  int M = F.size() + 1;
  int N = F[0].size() + 1;
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;

// Parallelizing the outer loop with OpenMP
#pragma omp parallel for num_threads(1)
  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < N - 1; j++) {
      // Node coordinates
      Point P_ij = {(i + 0.5) * h1, (j + 0.5) * h2};
      Point P_ij_diag = {P_ij.x + h1, P_ij.y + h2};

      F[i][j] = intersectionArea(P_ij, P_ij_diag);
    }
  }
}

double product(std::vector<std::vector<double>> &v1,
               std::vector<std::vector<double>> &v2, double h1, double h2) {
  double res = 0;
// Parallelizing the product calculation with OpenMP
#pragma omp parallel for reduction(+ : res) num_threads(1)
  for (int i = 0; i < static_cast<int>(v1.size()); i++) {
    for (int j = 0; j < static_cast<int>(v1[0].size()); j++) {
      res += h1 * h2 * v1[i][j] * v2[i][j];
    }
  }
  return res;
}

void calculateW(const std::vector<std::vector<double>> &a,
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

  for (int iter = 0; iter < maxIterations; iter++) {
    std::vector<std::vector<double>> newW = W;
    double maxChange = 0.0;

// Get residuals
#pragma omp parallel for num_threads(1)
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
#pragma omp parallel for num_threads(1)
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
#pragma omp parallel for num_threads(1)
    for (int i = 0; i < M - 1; i++) {
      for (int j = 0; j < N - 1; j++) {
        int I = i + 1;
        int J = j + 1;
        newW[I][J] = W[I][J] - theta * r[I][J];
// Track the maximum change
#pragma omp critical
        maxChange = std::max(
            maxChange, std::abs(newW[I][J] - W[I][J])); // fix: access to newW
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

void printABF(const std::vector<std::vector<double>> &a,
              const std::vector<std::vector<double>> &b,
              const std::vector<std::vector<double>> &F) {
  int M = F.size() + 1;
  int N = F[0].size() + 1;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N - 1; j++) {
      std::cout << " a[" << i + 1 << "][" << j + 1 << "] = " << a[i][j];
    }
    std::cout << std::endl;
  }

  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << " "
                << "b[" << i + 1 << "][" << j + 1 << "] = " << b[i][j];
    }
    std::cout << std::endl;
  }

  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < N - 1; j++) {
      std::cout << " F[" << i + 1 << "][" << j + 1 << "] = " << F[i][j];
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

int main() {
  int M = 10;
  int N = 10;

  // Example matrices (a, b, F) initialized with some values
  std::vector<std::vector<double>> a(M, std::vector<double>(N - 1, 0.0));
  std::vector<std::vector<double>> b(M - 1, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> F(M - 1, std::vector<double>(N - 1, 0.0));

  computeA(a);
  computeB(b);
  computeF(F);

  std::vector<std::vector<double>> W(M + 1, std::vector<double>(N + 1, 0.0));
  int maxIterations = 1000;
  double tolerance = 1e-5;
  calculateW(a, b, F, W, maxIterations, tolerance);

  std::cout << "Solution W:\n";
  // Output the result
  for (int i = 0; i <= M; i++) {
    for (int j = 0; j <= N; j++) {
      std::cout << std::fixed << std::setprecision(4) << W[i][j] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}