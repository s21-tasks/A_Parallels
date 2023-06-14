#pragma once

#include "matrix.h"
#include "../sub/thread_pool/thread_pool.h"
#include <vector>

namespace s21 {

  class SLAEGauss {
  public:

    SLAEGauss();
    SLAEGauss(unsigned threads);
    SLAEGauss(const Matrix<double>& m, const std::vector<double>& consts);
    SLAEGauss(const Matrix<double>& m, const std::vector<double>& consts, unsigned threads);

    SLAEGauss(SLAEGauss &gauss) = delete;
    ~SLAEGauss() = default;

    void set_equations(const Matrix<double>& m, const std::vector<double>& consts);

    std::vector<double> UsualExecute();
    std::vector<double> ParallelExecute();



  private:
    Matrix<double> matrix;
    std::vector<double> constants;
    std::vector<double> solution;
    int n = 0;
    ThreadPool pool;
  };

} // s21
