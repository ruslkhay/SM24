#pragma once
#include "inc.h"

/// @brief Structure of grid node
struct Point {
  double x, y;
};

/// @brief Check if the grid node belongs to trapezoid with nodes in
/// A(0,0), B(3,0), C(2,3), D(0,3).
/// @param p node to check
/// @return `True` if the node is in the trapezoid, `False` otherwise
bool isInTrapezoid(Point p);

/// @brief Check if the grid node belongs to enclosing rectangle with nodes in
/// A1(0,0), B1(3,0), A2(0,3), B2(0,3).
/// @param p node to check
/// @return `True` if the node is in the rectangle, `False` otherwise
bool isInRectangle(Point p);

/// @brief Calculate other coordinate of the point, laying on the side BC
/// of the trapezoid.
/// BC lays on y = -3x + 9 line
/// @param axis Either takes `x` or `y`
/// @param val Value of the coordinate, that corresponds to `axis` parameter
/// @return Second coordinate of the point
double lineBC(char axis, double val);

double horizontalShiftLen(Point left, Point right);

double verticalShiftLen(Point bottom, Point top);

/// @brief Calculate grid cell intersection area with trapezoid
/// @param bottomLeftPoint node in the bottom left corner of cell
/// @param topRightPoint node in the top right corner of cell
/// @return Area S_ij value
double intersectionArea(Point bottomLeftPoint, Point topRightPoint);

/// @brief Print results for verification
void printABF(const std::vector<std::vector<double>> &a,
              const std::vector<std::vector<double>> &b,
              const std::vector<std::vector<double>> &F);
