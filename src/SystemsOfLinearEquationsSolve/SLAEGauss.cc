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
    n = (int)constants.size();
    solution = std::vector<double>(n);
  }

  std::vector<double> SLAEGauss::UsualExecute() {
    for (int k = 0; k < n - 1; ++k) {
      StraightStrokeStep(k);
    }

    for (int i = n-1; i >=0;--i) {
      ReverseStrokeStep(i);
    }
    return solution;
  }

  std::vector<double> SLAEGauss::ParallelExecute() {
    return {};
  }

  void SLAEGauss::StraightStrokeStep(int k) {
    for (int i = k + 1; i < n; ++i) {
      double factor = matrix[i][k] / matrix[k][k];
      constants[i] -= factor * constants[k];
      for (int j = k; j < n; ++j) {
        matrix[i][j] -= factor * matrix[k][j];
      }
    }
  }

  void SLAEGauss::ReverseStrokeStep(int i) {
    double sum = 0;
    for (int j = i + 1; j < n; ++j)
      sum += matrix[i][j] * solution[j];
    solution[i] = (constants[i] - sum) / matrix[i][i];
  }

} // s21
