#pragma once

#include <vector>

namespace s21 {

/**
 * @struct TsmResult
 * @brief Represents the result of solving the Traveling Salesman Problem (TSP).
 * The TsmResult struct contains the information about the shortest path and its
 * length obtained from solving the Traveling Salesman Problem (TSP). It
 * includes a vector of vertices representing the order in which the vertices
 * should be visited to form the shortest path, and the distance representing
 * the total length of the path.
 */
struct TsmResult {
  std::vector<int> vertices;
  double distance = 0;

  /**
   * @brief Default constructor for TsmResult.
   */
  TsmResult() = default;

  /**
   * @brief Constructor for TsmResult.
   * @param d The total length of the shortest path.
   * @param v The vertices in the order they should be visited to form the
   * shortest path.
   */
  TsmResult(double d, const std::vector<int>& v) : distance(d), vertices(v) {}

  /**
   * @brief Assignment operator for TsmResult.
   * @param other The TsmResult object to assign from.
   * @return A reference to the assigned TsmResult object.
   */
  TsmResult& operator=(const TsmResult& other) {
    vertices = other.vertices;
    distance = other.distance;
    return *this;
  }

  /**
   * @brief Less-than comparison operator for TsmResult.
   * @param other The TsmResult object to compare against.
   * @return True if the distance of this TsmResult is less than the distance of
   * the other TsmResult, false otherwise.
   */
  bool operator<(const TsmResult& other) const {
    return distance < other.distance;
  }
};

}  // namespace s21