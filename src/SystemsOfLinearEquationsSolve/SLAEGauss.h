#pragma once

#include "matrix.h"
#include "../sub/thread_pool/thread_pool.h"
#include <vector>

namespace s21 {

  class SLAEGauss {
  public:

    SLAEGauss() = delete;
    SLAEGauss(const Matrix<double>& m, const std::vector<double>& consts);

    ~SLAEGauss() = default;

    void set_equations(const Matrix<double>& m, const std::vector<double>& consts);

    std::vector<double> UsualExecute();
    std::vector<double> ParallelExecute();



  private:

    void StraightStrokeStep(int k);
    void ReverseStrokeStep(int i);

    Matrix<double> matrix;
    std::vector<double> constants;
    std::vector<double> solution;
    int n = 0;




  };

} // s21
