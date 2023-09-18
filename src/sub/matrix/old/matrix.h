#pragma once

#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <vector>

namespace s21 {

template <class T>
class Matrix {
 private:
  std::unique_ptr<T[]> matrix_;
  size_t rows_, cols_;

  class RowProxy {
    Matrix<T> &matrix_;
    size_t row_;

   public:
    RowProxy(Matrix<T> &matrix, size_t row) : matrix_(matrix), row_(row) {}
    T &operator[](size_t col) {
      matrix_.row_col_error(row_, col, "operator[][]");
      return matrix_(row_, col);
    }
  };

  void error(const std::string &name) const {
    throw std::invalid_argument("Matrix: " + name + " error");
  }

  void row_col_error(size_t row, size_t col, std::string name) const {
    if (row >= rows_ || row < 0 || col >= cols_ || col < 0) error(name);
  }

 public:
  void Loop(const std::function<void(int, int)> &func) {
    for (size_t k = 0; k < rows_; ++k) {
      for (unsigned long g = 0; g < cols_; ++g) {
        func(k, g);
      }
    }
  }
  void Loop(const std::function<void(int, int)> &func) const {
    for (int k = 0; k < rows_; ++k) {
      for (int g = 0; g < cols_; ++g) {
        func(k, g);
      }
    }
  }

  Matrix() noexcept : rows_(0), cols_(0) {}

  Matrix(int rows, int cols)
      : rows_(rows), cols_(cols), matrix_(new T[rows * cols]()) {}

  Matrix(int square)
      : rows_(square), cols_(square), matrix_(new T[square * square]()) {}

  Matrix(int rows, int cols, const T &value)
      : rows_(rows), cols_(cols), matrix_(new T[rows * cols]) {
    Loop([&](int k, int g) { matrix_[k * cols_ + g] = value; });
  }

  Matrix(int rows, int cols, std::function<const T(void)> value_func)
      : matrix_(new T[rows * cols]), rows_(rows), cols_(cols) {
    Loop([&](int k, int g) { matrix_[k * cols_ + g] = value_func(); });
  }

  Matrix(int rows, int cols, std::function<void(T &)> value_func)
      : rows_(rows), cols_(cols), matrix_(new T[rows * cols]) {
    Loop([&](int k, int g) { value_func(matrix_[k * cols_ + g]); });
  }

  Matrix(int rows, int cols, std::function<void(T &, int, int)> value_func)
      : rows_(rows), cols_(cols), matrix_(new T[rows * cols]) {
    Loop([&](int k, int g) { value_func(matrix_[k * cols_ + g], k, g); });
  }

  Matrix(const Matrix &other)
      : matrix_(new T[other.rows_ * other.cols_]),
        rows_(other.rows_),
        cols_(other.cols_) {
    Loop([&](int k, int g) { matrix_[k * cols_ + g] = other(k, g); });
  }

  Matrix(std::initializer_list<std::initializer_list<T>> const &items)
      : rows_(items.size()),
        cols_(items.begin()->size()),
        matrix_(new T[items.size() * items.begin()->size()]) {
    int k = 0;
    for (auto const &row : items) {
      int g = 0;
      for (auto const &cell : row) {
        matrix_[k * cols_ + (g++)] = cell;
      }
      ++k;
    }
  }

  Matrix(Matrix &&other) noexcept
      : rows_(other.rows_),
        cols_(other.cols_),
        matrix_(std::move(other.matrix_)) {
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }

  RowProxy operator[](size_t row) { return RowProxy(*this, row); }

  T &operator()(size_t row, size_t col) {
    row_col_error(row, col, "operator()");
    return matrix_[row * cols_ + col];
  }

  const T &operator()(size_t row, size_t col) const {
    row_col_error(row, col, "const operator()");
    return matrix_[row * cols_ + col];
  }

  void Resize(size_t new_rows, size_t new_cols) {
    if (new_rows <= 0 || new_cols <= 0) error("Resize");
    Matrix<T> new_matrix(new_rows, new_cols);
    int min_rows = std::min(rows_, new_rows);
    int min_cols = std::min(cols_, new_cols);
    for (int k = 0; k < min_rows; ++k) {
      for (int g = 0; g < min_cols; ++g) {
        new_matrix(k, g) = matrix_[k * cols_ + g];
      }
    }
    *this = std::move(new_matrix);
  }

  void ResizeRows(size_t new_rows) { Resize(new_rows, cols_); }

  void ResizeCols(size_t new_cols) { Resize(rows_, new_cols); }

  Matrix &operator=(const Matrix &other) {
    if (this != &other) {
      Matrix temp(other);
      std::swap(rows_, temp.rows_);
      std::swap(cols_, temp.cols_);
      std::swap(matrix_, temp.matrix_);
    }
    return *this;
  }

  Matrix &operator=(Matrix &&other) noexcept {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
    return *this;
  }

  void SwapRows(int row1, int row2) {
    if (row1 == row2) return;
    if (row1 >= rows_ || row1 < 0 || row2 >= rows_ || row2 < 0)
      error("SwapRows");
    for (int col = 0; col < cols_; ++col) {
      std::swap(matrix_[row1 * cols_ + col], matrix_[row2 * cols_ + col]);
    }
  }

  void SwapCols(int col1, int col2) {
    if (col1 >= cols_ || col1 < 0 || col2 >= cols_ || col2 < 0)
      error("SwapCols");
    if (col1 == col2) return;
    for (int row = 0; row < rows_; ++row) {
      std::swap(matrix_[row * cols_ + col1], matrix_[row * cols_ + col2]);
    }
  }

  void ForEach(std::function<void(const T &)> func) const {
    Loop([&](int k, int g) { func(matrix_[k * cols_ + g]); });
  }

  void ForEach(std::function<void(const T &)> func) {
    Loop([&](int k, int g) { func(matrix_[k * cols_ + g]); });
  }

  void Fill(std::function<void(T &)> func) {
    Loop([&](int k, int g) { func(matrix_[k * cols_ + g]); });
  }

  const T &Get(int row, int col) const {
    row_col_error(row, col, "Get");
    return matrix_[row * cols_ + col];
  }

  void Set(int row, int col, T value) {
    row_col_error(row, col, "Set");
    matrix_[row * cols_ + col] = value;
  }

  size_t GetRows() const { return rows_; }

  size_t GetCols() const { return cols_; }

  template <class Mul>
  Matrix operator*(const Mul num) const {
    return Matrix(rows_, cols_, [&](T &cell, int k, int g) {
      cell = matrix_[k * cols_ + g] * num;
    });
  }

  template <class Mul>
  void operator*=(const Mul num) {
    Loop([&](int k, int g) { matrix_[k * cols_ + g] *= num; });
  }

  Matrix operator*(Matrix &other) const {
    if (cols_ != other.rows_) error("Mul");
    return Matrix(rows_, other.cols_, [&](T &cell, int k, int g) {
      cell = 0;
      for (int i = 0; i < cols_; ++i) {
        cell += matrix_[k * cols_ + i] * other(i, g);
      }
    });
  }

  void operator*=(Matrix &other) {
    auto temp = *this;
    *this = temp * other;
  }

  Matrix operator+(Matrix &other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) error("Sum");
    return Matrix(rows_, cols_, [&](T &cell, int k, int g) {
      cell = matrix_[k * cols_ + g] + other(k, g);
    });
  }

  void operator+=(Matrix &other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) error("Sum");
    Loop([&](int k, int g) { matrix_[k * cols_ + g] += other(k, g); });
  }

  Matrix operator-(Matrix &other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) error("Sub");
    return Matrix(rows_, cols_, [&](T &cell, int k, int g) {
      cell = matrix_[k * cols_ + g] - other(k, g);
    });
  }

  void operator-=(Matrix &other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) error("Sub");
    Loop([&](int k, int g) { matrix_[k * cols_ + g] -= other(k, g); });
  }

  Matrix Transpose() const {
    return Matrix(cols_, rows_, [&](T &cell, int k, int g) {
      cell = matrix_[g * cols_ + k];
    });
  }

  void TransposeThis() {
    auto temp = *this;
    *this = temp.Transpose();
  }

  bool operator==(const Matrix &other) const {
    if (cols_ != other.cols_ || rows_ != other.rows_) return false;
    for (int k = 0; k < rows_; ++k) {
      for (int g = 0; g < cols_; ++g) {
        if (other(k, g) != matrix_[k * cols_ + g]) {
          return false;
        }
      }
    }
    return true;
  }
  bool operator!=(Matrix &other) { return !operator==(other); }

  Matrix Minor(int row, int col) {
    row_col_error(row, col, "Minor");
    Matrix res(rows_ - 1, cols_ - 1);
    int i = 0, j;
    for (int k = 0; k < rows_; k++) {
      j = 0;
      if (k == row) continue;
      for (int g = 0; g < cols_; g++) {
        if (g != col) {
          res.matrix_[i * cols_ + j] = matrix_[k * cols_ + g];
          j++;
        }
      }
      i++;
    }
    return res;
  }
  T Determinant() {
    if (rows_ != cols_) error("Determinant");
    std::vector<std::vector<T>> lu(rows_, std::vector<T>(cols_));
    std::vector<int> pivot(rows_);
    Loop([&](int i, int j) { lu[i][j] = matrix_[i * cols_ + j]; });
    int sign = 1;
    for (int i = 0; i < rows_; i++) {
      pivot[i] = i;
    }
    for (int k = 0; k < rows_ - 1; k++) {
      T max_value = 0;
      int max_index = k;
      for (int i = k; i < rows_; i++) {
        T abs_value = std::abs(lu[i][k]);
        if (abs_value > max_value) {
          max_value = abs_value;
          max_index = i;
        }
      }
      if (max_index != k) {
        std::swap(lu[k], lu[max_index]);
        std::swap(pivot[k], pivot[max_index]);
        sign *= -1;
      }
      for (int i = k + 1; i < rows_; i++) {
        T factor = lu[i][k] / lu[k][k];
        for (int j = k + 1; j < cols_; j++) {
          lu[i][j] -= factor * lu[k][j];
        }
        lu[i][k] = factor;
      }
    }
    T det = sign;
    for (int i = 0; i < rows_; i++) {
      det *= lu[i][i];
    }
    return det;
  }

  friend std::ostream &operator<<(std::ostream &os, const Matrix<T> &m) {
    for (int i = 0; i < m.rows_; ++i) {
      for (int j = 0; j < m.cols_; ++j) {
        os << m(i, j) << '\t';
      }
      os << "\n";
    }
    return os;
  }
  void Print() const { std::cout << *this; }

  virtual ~Matrix() {}
};

}  // namespace s21
