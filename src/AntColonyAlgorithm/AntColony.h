#pragma once

#include <limits>

#include "Helpers.h"
#include "../sub/Matrix/Matrix.h"
#include "../sub/Utility/Utility.h"

namespace s21 {

template <class T>
class Ant;

/**
 * @class AntColony
 * @brief Implements the Ant Colony Optimization algorithm.
 *
 * The AntColony class represents an Ant Colony Optimization algorithm for
 * solving optimization problems on a graph represented by a matrix. It uses a
 * colony of ants to iteratively construct solutions and update pheromone levels
 * on graph edges. The class provides methods for initializing and running the
 * algorithm, as well as accessing the result.
 * @tparam T The type of the matrix elements.
 */
template <class T>
class AntColony {
 public:
  /**
   * @brief Constructor for the AntColony class.
   *
   * Constructs an AntColony object with the specified parameters.
   *
   * @param graph The graph represented as a matrix.
   * @param alpha The alpha parameter, pheromone influence (default: 0.8).
   * @param beta The beta parameter, edge weight influence (default: 1.2).
   * @param rho The rho parameter, intensity of evaporation of pheromones
   * (default: 0.5).
   * @param iterations The number of iterations to perform (default: 500).
   * @param ants_count_k The number of ants in the colony is product of
   * ants_count_k and graph size (default: 1).
   * @param Q_k The Q parameter, the number of pheromones in one ant is product
   * of the Q_k to the average weight of the edge (default: 0.6).
   * @param default_pheromone The default pheromone level on edges (default:
   * 0.6).
   */
  AntColony(const Matrix<T> &graph, double alpha = 0.8, double beta = 1.2,
            double rho = 0.5, int iterations = 500, int ants_count_k = 1,
            double Q_k = 0.6, double defult_pheromone = 0.6);

  /**
   * @brief Solve the optimization problem using the Ant Colony Optimization
   * algorithm.
   *
   * This method runs the Ant Colony Optimization algorithm for the specified
   * number of iterations and returns the result as a TsmResult object, which
   * contains the best route found by the algorithm.
   *
   * @return The result of the Ant Colony Optimization algorithm as a TsmResult
   * object.
   */
  TsmResult Solve();

 protected:
  friend class Ant<T>;
  double alpha_, beta_, rho_, Q_;
  int iterations_, size_;
  const Matrix<T> &graph_;
  Matrix<double> pheromone_;
  std::vector<Ant<T>> ants_;
  TsmResult result_;

  void LeavePheromone(const TsmResult &route);
  void CreateAnts();
  virtual void AntsRun();
  void PheromoneEvaporation();
  void PheromoneUpdate();
};

template <class T>
class Ant {
 public:
  Ant(const AntColony<T> &colony);
  Ant(const Ant &other) : Ant(other.colony_) {}

  void Run();
  void Clear();

  const TsmResult &GetRoute() const {
    return route_;
  }

 private:
  friend class AntColony<T>;
  std::vector<double> probabilities_;
  std::vector<bool> memory_;
  int position_, visited_ = 0, start_;
  TsmResult route_;
  const AntColony<T> &colony_;

  void Visit(int move);
  double CountProbabilities();
  bool RandomVisit(double sum);
};

template <class T>
AntColony<T>::AntColony(const Matrix<T> &graph, double alpha, double beta,
                        double rho, int iterations, int ants_count_k,
                        double Q_k, double defult_pheromone)
    : graph_(graph),
      alpha_(alpha),
      beta_(beta),
      rho_(1.0 - rho),
      size_(graph.GetCols()),
      ants_(ants_count_k * size_, Ant<T>(*this)),
      iterations_(iterations),
      Q_(0.0),
      pheromone_(size_, size_, defult_pheromone),
      result_(INFINITY, std::vector<int>(size_ + 1, 0)) {
  
  // for (int k = 0; k < ants_count_k * size_; ++k) {
  //   ants_.emplace_back(*this);
  // }

  graph_.ForEach([&](const T &cell) { Q_ += cell; });
  Q_ *= Q_k / (size_ * size_);
}

template <class T>
TsmResult AntColony<T>::Solve() {
  CreateAnts();
  while (--iterations_ != 0) {
    AntsRun();
    PheromoneEvaporation();
    PheromoneUpdate();
  }
  return result_;
}

template <class T>
void AntColony<T>::CreateAnts() {
  for (int k = 0; k < ants_.size(); ++k) {
    ants_[k].start_ = k % size_;
  }
}

template <class T>
void AntColony<T>::AntsRun() {
  for (auto &ant : ants_) {
    ant.Run();
    if (ant.route_.distance < result_.distance) {
      result_ = ant.route_;
    }
  }
}

template <class T>
void AntColony<T>::PheromoneEvaporation() {
  pheromone_.Fill([&](double &cell) { cell *= rho_; });
}

template <class T>
void AntColony<T>::PheromoneUpdate() {
  for (auto &ant : ants_) {
    LeavePheromone(ant.route_);
    ant.Clear();
  }
}

template <class T>
void AntColony<T>::LeavePheromone(const TsmResult &route) {
  auto iter = route.vertices.begin();
  int prev = *iter;
  double delta = Q_ / route.distance;
  while (++iter != route.vertices.end()) {
    pheromone_(prev, *iter) += delta;
    prev = *iter;
  }
}

template <class T>
Ant<T>::Ant(const AntColony<T> &colony)
    : memory_(colony.size_, false),
      colony_(colony),
      probabilities_(colony.size_, -1.0),
      route_(0, std::vector<int>(colony.size_ + 1, 0)) {}

template <class T>
void Ant<T>::Run() {
  position_ = start_;
  memory_[start_] = true;
  route_.vertices[0] = start_;

  while (visited_ < colony_.size_ - 1) {
    if (!RandomVisit(CountProbabilities())) {
      route_.distance = std::numeric_limits<double>::infinity();
      return;
    }
  }

  if (colony_.graph_.Get(route_.vertices[visited_], start_) == 0) {
    route_.distance = std::numeric_limits<double>::infinity();
  } else {
    route_.distance += colony_.graph_.Get(position_, start_);
    route_.vertices[visited_ + 1] = start_;
  }
}

template <class T>
double Ant<T>::CountProbabilities() {
  double sum = 0.0;
  for (int i = 0; i < memory_.size(); ++i) {
    if (colony_.graph_.Get(position_, i) != 0 && !memory_[i]) {
      probabilities_[i] =
          pow(colony_.pheromone_.Get(position_, i), colony_.alpha_) *
          pow(1.0 / colony_.graph_.Get(position_, i), colony_.beta_);
      sum += probabilities_[i];
    } else {
      probabilities_[i] = 0.0;
    }
  }
  return sum;
}

template <class T>
bool Ant<T>::RandomVisit(double sum) {
  double selection = Random::Uniform<double>(0, sum);
  for (int i = 0; i < probabilities_.size(); ++i) {
    selection -= probabilities_[i];
    if (selection <= 0.0 && probabilities_[i] > 0.0) {
      Visit(i);
      return true;
    }
  }
  return false;
}

template <class T>
void Ant<T>::Clear() {
  memory_.assign(memory_.size(), false);
  visited_ = 0;
  route_.distance = 0;
}

template <class T>
void Ant<T>::Visit(int move) {
  route_.distance += colony_.graph_.Get(position_, move);
  position_ = move;
  memory_[move] = true;
  route_.vertices[++visited_] = move;
}

}  // namespace s21
