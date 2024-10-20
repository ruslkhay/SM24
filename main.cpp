#include <cmath>
#include <iomanip> // for std::setprecision
#include <iostream>
#include <vector>

const double eps = 0.01;

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

// Function to compute the length of the segment in TR
double lengthInTrapezoid(Point p1, Point p2) {
  // if (isInRectangle(p1) && isInRectangle(p2))
  bool horizShift = p1.y == p2.y;
  bool vertShift = p1.x == p2.x;
  if (isInTrapezoid(p1) && isInTrapezoid(p2)) {
    if (vertShift) {
      return std::abs(p2.y - p1.y);
    }
    if (horizShift) {
      return std::abs(p2.x - p1.x);
    }
  } else if (isInTrapezoid(p1) || isInTrapezoid(p2)) {
    // If only one point is inside the trapezoid, calculate the distance to the
    // edge
    Point inPoint = isInTrapezoid(p1) ? p1 : p2;
    Point outPoint = isInTrapezoid(p1) ? p2 : p1;

    if (!isInRectangle(outPoint)) {
      if (vertShift) {
        // If point is in sub-rectangular of trapezoid
        if (inPoint.x <= 2) {
          return (outPoint.y - inPoint.y) * 1 / 2;
        }
        // If point is in sub-triangle of trapezoid
        else {
          double edgeY = -3 * outPoint.x + 9;
          return edgeY - inPoint.y;
        }
      }
      if (horizShift) {
        double edgeX = (9 - outPoint.y) / 3;
        return edgeX - inPoint.x;
      }
    } else {
      if (vertShift) {
        double edgeY = -3 * outPoint.x + 9;
        return edgeY - inPoint.y;
      }
      if (horizShift) {
        double edgeX = (9 - outPoint.y) / 3;
        return edgeX - inPoint.x;
      }
    }
  }
  return 0; // If both points are outside and not in trapezoid
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
      p1 = {(9 - topRightPoint.y) / 3, topLeftPoint.y};
    } else {
      p1 = {topLeftPoint.x, -3 * topLeftPoint.x + 9};
    }
    if (isInTrapezoid(bottomRightPoint)) {
      p2 = {bottomRightPoint.x, -3 * bottomRightPoint.x + 9};
    } else {
      p2 = {(9 - bottomRightPoint.y) / 3, bottomRightPoint.y};
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
      std::cout << bottomLeftPoint.x << ',' << bottomLeftPoint.y << ' '
                << topRightPoint.x << ',' << topRightPoint.y << " is a triangle"
                << std::endl;
      auto sideV = p1.y - bottomLeftPoint.y;
      auto sideH = p2.x - bottomLeftPoint.x;
      S_ij = sideV * sideH * 1 / 2;
    }
    // If gird rectangle without interception (complement of interception) is a
    // triangle
    if ((p1OnHor | p1InCorner) && (p2OnVert | p2InCorner)) {
      std::cout << bottomLeftPoint.x << ',' << bottomLeftPoint.y << ' '
                << topRightPoint.x << ',' << topRightPoint.y
                << " complement is a triangle" << std::endl;
      auto sideV = topRightPoint.x - p1.x;
      auto sideH = topRightPoint.y - p2.y;
      S_ij = h1 * h2 - sideV * sideH * 1 / 2;
    }
    // If S_ij is a trapezoid with a vertical bases
    if (p1OnVert && p2OnVert) {
      std::cout << bottomLeftPoint.x << ',' << bottomLeftPoint.y << ' '
                << topRightPoint.x << ',' << topRightPoint.y
                << " is a trapezoid with vertical bases" << std::endl;
      auto h = bottomRightPoint.x - bottomLeftPoint.x;
      auto base1 = p2.y - bottomRightPoint.y;
      auto base2 = p1.y - bottomLeftPoint.y;
      S_ij = h * (base1 + base2) / 2;
    }
    // If S_ij is a trapezoid with a horizontal bases
    if (p1OnHor && p2OnHor) {
      std::cout << bottomLeftPoint.x << ',' << bottomLeftPoint.y << ' '
                << topRightPoint.x << ',' << topRightPoint.y
                << " is a trapezoid with a horizontal bases" << std::endl;
      auto h = topLeftPoint.y - bottomLeftPoint.y;
      auto base1 = p1.x - topLeftPoint.x;
      auto base2 = p2.x - bottomLeftPoint.x;
      S_ij = h * (base1 + base2) / 2;
    }
  }
  return S_ij / (h1 * h2);
}

void computeGrid(std::vector<std::vector<double>> &a,
                 std::vector<std::vector<double>> &b,
                 std::vector<std::vector<double>> &F) {

  int M = a.size();
  int N = a[0].size();
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;

  // Calculate values for each inner node
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      // Node coordinates
      Point P_ij = {((i + 1) - 0.5) * h1, ((j + 1) - 0.5) * h2};
      Point P_ij_next = {P_ij.x, P_ij.y + h2};
      Point P_ji_next = {P_ij.x + h1, P_ij.y};

      double lv_ij = lengthInTrapezoid(P_ij, P_ij_next);
      double lh_ij = lengthInTrapezoid(P_ij, P_ji_next);

      if (lv_ij == h2) {
        a[i][j] = 1;
      } else {
        a[i][j] = 1 / h2 * lv_ij + (1 - lv_ij / h2) * (1 / eps);
      }

      if (lh_ij == h1) {
        b[i][j] = 1;
      } else {
        b[i][j] = 1 / h1 * lh_ij + (1 - lh_ij / h1) * (1 / eps);
      }

      if (i != M - 1 && j != N - 1) {
        Point P_ij_diag = {P_ij.x + h1, P_ij.y + h2};
        F[i][j] = intersectionArea(P_ij, P_ij_diag);
      }
    }
  }
}

void calculateW(const std::vector<std::vector<double>> &a,
                const std::vector<std::vector<double>> &b,
                const std::vector<std::vector<double>> &F,
                std::vector<std::vector<double>> &W, int maxIterations,
                double tolerance) {

  int M = a.size();
  int N = a[0].size();
  std::vector<std::vector<double>> r(M + 1, std::vector<double>(N + 1, 0.0));
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;
  // Perform the iterative steepest descent
  for (int iter = 0; iter < maxIterations; iter++) {
    std::vector<std::vector<double>> newW = W;
    double maxChange = 0.0;

    // Iterate through the elements to get residuals
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
    //  of W (excluding the boundary conditions)
    for (int i = 0; i < M - 1; i++) {
      for (int j = 0; j < N - 1; j++) {
        int I = i + 1;
        int J = j + 1;
        // Calculate the finite difference terms
        double term1 = (a[i][j] * (r[I][J] - r[I - 1][J])) / (h1 * h1);
        double term2 = (a[i + 1][j] * (r[I + 1][J] - r[I][J])) / (h1 * h1);
        double term3 = (b[i][j] * (r[I][J] - r[I][J - 1])) / (h2 * h2);
        double term4 = (b[i][j + 1] * (r[I][J + 1] - r[I][J])) / (h2 * h2);

        // Compute the new value of W based on the equation
        double theta = r[I][J] / (term1 - term2 + term3 - term4);
        newW[I][J] = W[I][J] - theta * r[I][J];
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

int main() {
  // Define trapezoid vertices
  // BC side lays on y = -3x + 9 line
  // const Point A = {0, 0};
  // const Point B = {3, 0};
  // const Point C = {2, 3};
  // const Point D = {0, 3};

  // Example usage
  int M = 5;
  int N = 5;

  // Example matrices (a, b, F) initialized with some values
  std::vector<std::vector<double>> a(M, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> b(M, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> F(M - 1, std::vector<double>(N - 1, 0.0));
  computeGrid(a, b, F);

  // Print results for verification
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << "a[" << i << "][" << j << "] = " << a[i][j] << ", "
                << "b[" << i << "][" << j << "] = " << b[i][j];
      if (i != M - 1 && j != N - 1) {
        std::cout << ", "
                  << "F[" << i << "][" << j << "] = " << F[i][j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  std::vector<std::vector<double>> W(M + 1, std::vector<double>(N + 1, 0.0));
  int maxIterations = 1000;
  double tolerance = 1e-5;
  calculateW(a, b, F, W, maxIterations, tolerance);

  // Output the result
  for (int i = 0; i <= M; i++) {
    for (int j = 0; j <= N; j++) {
      std::cout << std::fixed << std::setprecision(4) << W[i][j] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}