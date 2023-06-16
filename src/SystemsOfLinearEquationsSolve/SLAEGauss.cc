#include "SLAEGauss.h"

namespace s21 {

  SLAEGauss::SLAEGauss() = default;




  SLAEGauss::SLAEGauss(const Matrix<double> &m, const std::vector<double> &consts){
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
      for (int i = k + 1; i < n; ++i) {
        double factor = matrix[i][k] / matrix[k][k];
        constants[i] -= factor * constants[k];
        for (int j = k; j < n; ++j) {
          matrix[i][j] -= factor * matrix[k][j];
        }
      }
    }

    for (int i = n-1; i >=0;--i) {
      double sum = 0;
      for (int j = i + 1; j < n; ++j)
        sum += matrix[i][j] * solution[j];
      solution[i] = (constants[i] - sum) / matrix[i][i];
    }
    return solution;
  }


  std::vector<double> SLAEGauss::ParallelExecute(unsigned threads) {
    ThreadPool pool(threads);
    for (int k = 0; k < n - 1; ++k) {
      for (int i = k + 1; i < n; ++i) {
        double factor = matrix[i][k] / matrix[k][k];
        constants[i] -= factor * constants[k];
        for (int j = k; j < n; ++j) {
          pool.AddTask([=]() { matrix(i, j) -= factor * matrix(k, j); });
        }
        pool.WaitForComplete();
      }
    }

    for (int i = n - 1; i >= 0; --i) {
      double sum = 0;
      for (int j = i + 1; j < n; ++j)
        pool.AddTask([=, &sum]() { sum += matrix(i, j) * solution[j]; });
      pool.WaitForComplete();
      solution[i] = (constants[i] - sum) / matrix[i][i];
    }
    return solution;
  }


  void SLAEGauss::Print() {
    for (int i = 0;  i < matrix.GetRows(); ++i) {
      for (int j = 0;  j < matrix.GetCols(); ++j) {
        std::cout << (matrix[i][j] >= 0 ? " + " : " - ");
        std::cout << std::fabs(matrix[i][j]) << " * " << static_cast<char>(97+i) ;
      }
      std::cout << " = " << constants[i] << std::endl;
    }
  }


} // s21
