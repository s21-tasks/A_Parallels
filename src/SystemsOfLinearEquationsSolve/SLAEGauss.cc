//
// Created by Иван Захаров on 04.06.2023.
//

#include "SLAEGauss.h"

namespace s21 {

  SLAEGauss::SLAEGauss(const Matrix<double> &m, const std::vector<double> &consts) {
    set_equations(m,consts);
  }

  void SLAEGauss::set_equations(const Matrix<double>& m, const std::vector<double>& consts) {
    if (m.GetRows() != consts.size() || m.GetRows()!= m.GetCols())
      throw std::out_of_range("invalid sizes!");
    matrix = m;
    constants = consts;
  }

  std::vector<double> SLAEGauss::UsualExecute() {
    const int n = (int)constants.size();
    for (int k = 0; k < n - 1; ++k) {
      for (int i = k + 1; i < n; ++i) {
        double factor = matrix[i][k] / matrix[k][k];
        constants[i] -= factor * constants[k];
        for (int j = k; j < n; ++j) {
          matrix[i][j] -= factor * matrix[k][j];
        }
      }
    }

    std::vector<double> solution(n);

    for (int i = n-1; i >=0;--i) {
      double sum = 0;
      for (int j = i + 1; j < n; ++j)
        sum += matrix[i][j] * solution[j];
      solution[i] = (constants[i] - sum) / matrix[i][i];
    }
    return solution;
  }

  std::vector<double> SLAEGauss::ParallelExecute() {
    return {};
  }


} // s21