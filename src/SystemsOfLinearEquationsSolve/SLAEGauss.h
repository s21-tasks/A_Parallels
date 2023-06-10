#pragma once

#include "matrix.h"
#include <vector>

namespace s21 {

  class SLAEGauss {
  public:

    SLAEGauss() = delete;
    SLAEGauss(const Matrix<double>& m, const std::vector<double>& consts);


    void set_equations(const Matrix<double>& m, const std::vector<double>& consts);

    std::vector<double> UsualExecute();
    std::vector<double> ParallelExecute();



  private:
    Matrix<double> matrix;
    std::vector<double> constants;




  };

} // s21
