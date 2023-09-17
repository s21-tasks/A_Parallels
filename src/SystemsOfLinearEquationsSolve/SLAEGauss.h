#pragma once

#include <thread>
#include <vector>

#include "matrix.h"

namespace s21 {

class SLAEGauss {
 public:
  SLAEGauss();
  SLAEGauss(const Matrix<double>& m, const std::vector<double>& consts);

  SLAEGauss(SLAEGauss& gauss) = default;
  ~SLAEGauss() = default;

  void set_equations(const Matrix<double>& m,
                     const std::vector<double>& consts);

  std::vector<double> UsualExecute();
  std::vector<double> ParallelExecute(unsigned threads);

  void Print();

 private:
  Matrix<double> matrix;
  std::vector<double> constants;
  std::vector<double> solution;
  int n = 0;
};

}  // namespace s21
