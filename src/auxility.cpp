#include "auxility.hpp"

bool isInTrapezoid(Point p) {
  return (p.y <= -3 * p.x + 9 && p.x >= 2 && p.x <= 3) ||
         (p.y >= 0 && p.y <= 3 && p.x >= 0 && p.x <= 2);
}

bool isInRectangle(Point p) {
  return p.y >= 0 && p.y <= 3 && p.x >= 0 && p.x <= 3;
}

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
  } else if (isInTrapezoid(left)) {
    auto rLimit = lineBC('y', right.y);
    return rLimit - left.x;
  } else {
    return 0;
  }
}

double verticalShiftLen(Point bottom, Point top) {
  if (isInTrapezoid(top)) {
    return top.y - bottom.y;
  } else if (isInTrapezoid(bottom)) {
    auto tLimit = lineBC('x', top.x);
    return tLimit - bottom.y;
  } else {
    return 0;
  }
}

double intersectionArea(Point bottomLeftPoint, Point topRightPoint) {
  double h1 = topRightPoint.x - bottomLeftPoint.x;
  double h2 = topRightPoint.y - bottomLeftPoint.y;
  double S_ij = 0;
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
    // If gird rectangle without interception (complement of interception) is a
    // triangle
    if ((p1OnHor | p1InCorner) && (p2OnVert | p2InCorner)) {
      auto sideV = topRightPoint.x - p1.x;
      auto sideH = topRightPoint.y - p2.y;
      S_ij = h1 * h2 - sideV * sideH * 1 / 2;
    }
    // If S_ij is a trapezoid with a vertical bases
    if (p1OnVert && p2OnVert) {
      auto h = bottomRightPoint.x - bottomLeftPoint.x;
      auto base1 = p2.y - bottomRightPoint.y;
      auto base2 = p1.y - bottomLeftPoint.y;
      S_ij = h * (base1 + base2) / 2;
    }
    // If S_ij is a trapezoid with a horizontal bases
    if (p1OnHor && p2OnHor) {
      auto h = topLeftPoint.y - bottomLeftPoint.y;
      auto base1 = p1.x - topLeftPoint.x;
      auto base2 = p2.x - bottomLeftPoint.x;
      S_ij = h * (base1 + base2) / 2;
    }
  }
  return S_ij;
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