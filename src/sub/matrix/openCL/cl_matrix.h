#pragma once

#include <algorithm>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <vector>

#define CL_HPP_TARGET_OPENCL_VERSION 300
#include <clBLAS.h>

#include <CL/opencl.hpp>

namespace CL {

class Arithmetic;

template <class T>
class Matrix {
  friend class Arithmetic;
  size_t cols_, rows_;
  cl::Context context_;
  cl::Device default_device_;
  cl::Buffer matrix_;

 public:
  Matrix(size_t row, size_t col)
      : rows_(row),
        cols_(col),
        context_(CL_DEVICE_TYPE_GPU),
        default_device_(cl::Device::getDefault()),
        matrix_(
            cl::Buffer(context_, CL_MEM_READ_WRITE, sizeof(T) * row * col)) {}

  Matrix(size_t row, size_t col, const std::vector<T> &vec) : Matrix(row, col) {
    cl::CommandQueue queue(context_, default_device_);
    queue.enqueueWriteBuffer(matrix_, CL_TRUE, 0, sizeof(T) * row * col,
                             vec.data());
  }
  Matrix(size_t row, size_t col, std::function<T(void)> func)
      : Matrix(row, col) {
    std::vector<T> vec(row * col);
    std::generate(vec.begin(), vec.end(), func);
    cl::CommandQueue queue(context_, default_device_);
    queue.enqueueWriteBuffer(matrix_, CL_TRUE, 0, sizeof(T) * row * col,
                             vec.data());
  }

  Matrix(std::initializer_list<std::initializer_list<T>> const &items)
      : Matrix(items.size(), items.begin()->size()) {
    cl::CommandQueue queue(context_, default_device_);
    size_t k = 0;
    std::vector<T> buffer(rows_ * cols_);
    for (auto const &row : items) {
      for (auto const &cell : row) {
        buffer[k++] = cell;
      }
    }
    queue.enqueueWriteBuffer(matrix_, CL_TRUE, 0, sizeof(T) * rows_ * cols_,
                             buffer.data());
  }

  void Update(const std::vector<T> &vec) {
    cl::CommandQueue queue(context_, default_device_);
    queue.enqueueWriteBuffer(matrix_, CL_TRUE, 0, sizeof(T) * rows_ * cols_,
                             vec.data());
  }

  const size_t GetCols() const { return cols_; }
  const size_t GetRows() const { return rows_; }

  void Set(size_t row, size_t col, T value) {
    cl::CommandQueue queue(context_, default_device_);
    queue.enqueueWriteBuffer(matrix_, CL_TRUE, sizeof(T) * (row * cols_ + col),
                             sizeof(T), &value);
  }
  T Get(size_t row, size_t col) const {
    cl::CommandQueue queue(context_, default_device_);
    T element;
    queue.enqueueReadBuffer(matrix_, CL_TRUE, sizeof(T) * (row * cols_ + col),
                            sizeof(T), &element);
    return element;
  }

  std::vector<T> ToVector() {
    std::vector<T> vec(rows_ * cols_);
    cl::CommandQueue queue(context_, default_device_);
    queue.enqueueReadBuffer(matrix_, CL_TRUE, 0, sizeof(T) * vec.size(),
                            vec.data());
    return vec;
  }
  const std::vector<T> ToVector() const {
    std::vector<T> vec(rows_ * cols_);
    cl::CommandQueue queue(context_, default_device_);
    queue.enqueueReadBuffer(matrix_, CL_TRUE, 0, sizeof(T) * vec.size(),
                            vec.data());
    return vec;
  }

  void Print() {
    auto vec = ToVector();
    for (int k = 0; k < rows_; ++k) {
      for (int g = 0; g < cols_; ++g) {
        std::cout << vec[cols_ * k + g] << ' ';
      }
      std::cout << '\n';
    }
  }
};

class Arithmetic {
 public:
  template <class T>
  static void Mul(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C,
                  const T alpha = 1.0, const double beta = 0.0) {
    cl::CommandQueue queue(A.context_, A.default_device_);
    cl_command_queue raw_queue = queue.get();
    // clblasSetup();
    functions<T>::mul(clblasRowMajor, clblasNoTrans, clblasNoTrans, A.rows_,
                      B.cols_, A.cols_, alpha, A.matrix_(), 0, A.cols_,
                      B.matrix_(), 0, B.cols_, beta, C.matrix_(), 0, C.cols_, 1,
                      &raw_queue, 0, nullptr, nullptr);
    // clblasTeardown();
    queue.finish();
  }
  template <class T>
  static void Sigmoid(Matrix<T> &destination, const Matrix<T> &biases) {}

 private:
  template <class T>
  struct functions;

  // static void clblasSetup() {
  //     static std::once_flag flag;
  //     std::call_once(flag, []() { clblasSetup(); });
  // }

  // static void clblasTeardown() {
  //     static std::once_flag flag;
  //     std::call_once(flag, []() { clblasTeardown(); });
  // }
};

template <>
struct Arithmetic::functions<float> {
  constexpr static auto mul = clblasSgemm;
};

template <>
struct Arithmetic::functions<double> {
  constexpr static auto mul = clblasDgemm;
};

}  // namespace CL
