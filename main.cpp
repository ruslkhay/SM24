#include <cmath>
#include <iostream>
#include <vector>

const double eps = 0.01;

struct Point {
  double x, y;
};

// Define trapezoid vertices
// BC side lays on y = -3x + 9 line
const Point A = {0, 0};
const Point B = {3, 0};
const Point C = {2, 3};
const Point D = {0, 3};

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
/// diagonal points
/// @param bottomLeftPoint bottom left corner
/// @param topRightPoint top right corner
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
      S_ij = (base1 + base2) / h;
    }
    // If S_ij is a trapezoid with a horizontal bases
    if (p1OnHor && p2OnHor) {
      std::cout << bottomLeftPoint.x << ',' << bottomLeftPoint.y << ' '
                << topRightPoint.x << ',' << topRightPoint.y
                << " is a trapezoid with a horizontal bases" << std::endl;
      auto h = topLeftPoint.y - bottomLeftPoint.y;
      auto base1 = p1.x - topLeftPoint.x;
      auto base2 = p2.y - bottomLeftPoint.y;
      S_ij = (base1 + base2) / h;
    }
  }
  return S_ij / (h1 * h2);
}

void computeGrid(int M, int N) {
  double h1 = 3.0 / M;
  double h2 = 3.0 / N;

  std::vector<std::vector<double>> a(M + 1, std::vector<double>(N + 1, 0));
  std::vector<std::vector<double>> b(M + 1, std::vector<double>(N + 1, 0));
  std::vector<std::vector<double>> F(M, std::vector<double>(N, 0));

  // Calculate values for each inner node
  for (int i = 1; i <= M; i++) {
    for (int j = 1; j <= N; j++) {
      // Node coordinates
      Point P_ij = {(i - 0.5) * h1, (j - 0.5) * h2};
      Point P_ij_next = {P_ij.x, P_ij.y + h2};
      Point P_ji_next = {P_ij.x + h1, P_ij.y};

      double lv_ij = lengthInTrapezoid(P_ij, P_ij_next);
      double lh_ij = lengthInTrapezoid(P_ij, P_ji_next);

      if (lv_ij == h2) {
        a[i][j] = 1;
      } else {
        a[i][j] = 1 / h2 * lv_ij + (1 - lv_ij / h2) * (1 / eps);
        // std::cout << "i=" << i << ' ' << "j=" << j << ' '
        //     << "lv=" << lv_ij << ' ' << "lh=" << lh_ij << ' '
        //     << "P_ij=" << P_ij.x << ',' << P_ij.y << std::endl;
      }

      if (lh_ij == h1) {
        b[i][j] = 1;
      } else {
        b[i][j] = 1 / h1 * lh_ij + (1 - lh_ij / h1) * (1 / eps);
      }

      if (i != M && j != N) {
        Point P_ij_diag = {P_ij.x + h1, P_ij.y + h2};
        F[i][j] = intersectionArea(P_ij, P_ij_diag);
      }
    }
  }

  // Print results for verification
  for (int i = 1; i <= M; i++) {
    for (int j = 1; j <= N; j++) {
      std::cout << "a[" << i << "][" << j << "] = " << a[i][j] << ", "
                << "b[" << i << "][" << j << "] = " << b[i][j];
      if (i != M && j != N) {
        std::cout << ", "
                  << "F[" << i << "][" << j << "] = " << F[i][j];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

int main() {
  int M = 6, N = 6; // Example input
  computeGrid(M, N);

  return 0;
}