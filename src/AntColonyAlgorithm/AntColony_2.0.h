#pragma once

#include "../thread/m_thread.h"
#include "AntColony.h"

namespace s21 {

template <class T>
class AntColony20 : public AntColony<T> {
  using AC = AntColony<T>;

 public:
  AntColony20(const Matrix<T> &graph, int threads = 4, double alpha = 0.8,
              double beta = 1.2, double rho = 0.5, int iterations = 500,
              int ants_count_k = 1, double Q_k = 0.6,
              double defult_pheromone = 0.6)
      : AC(graph, alpha, beta, rho, iterations, ants_count_k, Q_k,
           defult_pheromone),
        thm_(threads) {
    thm_.SetLoopFunction(AntColony<T>::ants_.size(),
                         [&](int k) { AC::ants_[k].Run(); });
  }

 private:
  int count = 0;

  ThreadManager thm_;

  void AntsRun() override {
    // thm_.Run();

    // thm_.LoopExecute(AntColony<T>::ants_.size(), [&] (int k) {
    //     AC::ants_[k].Run();
    // });

    if (++count == 1) {
      thm_.SetLoopFunction(AntColony<T>::ants_.size(), [&](int k) {
        std::cout << "FUNC\n";
        AC::ants_[k].Run();
      });
    }
    // thm_.SetLoopFunction(AntColony<T>::ants_.size(),
    // std::bind(&AntColony20<T>::ThreadFunc, this, std::placeholders::_1));
    thm_.Run();

    for (const auto &ant : AC::ants_) {
      if (ant.GetRoute().distance < AC::result_.distance) {
        AC::result_ = ant.GetRoute();
      }
    }
  }

  void ThreadFunc(int k) { AC::ants_[k].Run(); }
};

}  // namespace s21
