#pragma once

#include "matrix.h"
#include <cblas.h>

namespace s21 {

template<class T>
class BLAS {

    struct func;

    public:
        // C = alpha*A*B + beta*C
        static void Mul(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C, const T alpha = 1.0, const double beta = 0.0) {
            func::mul(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.rows_, B.cols_, A.cols_,
                alpha, A.Data(), A.cols_, B.Data(), B.cols_, beta, C.Data(), B.cols_);
        }
        static void MulAT(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C, const T alpha = 1.0, const double beta = 0.0) {
            func::mul(CblasRowMajor, CblasTrans, CblasNoTrans, A.cols_, B.cols_, A.rows_,
                alpha, A.Data(), A.cols_, B.Data(), B.cols_, beta, C.Data(), B.cols_);
        }
        static void MulBT(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C, const T alpha = 1.0, const double beta = 0.0) {
            func::mul(CblasRowMajor, CblasNoTrans, CblasTrans, A.rows_, B.rows_, A.cols_,
                alpha, A.Data(), A.cols_, B.Data(), B.cols_, beta, C.Data(), B.rows_);
        }
        static void Mul(const std::vector<T> &A, const Matrix<T> &B, Matrix<T> &C, const T alpha = 1.0, const double beta = 0.0) {
            func::mul(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, B.cols_, B.rows_,
                alpha, A.data(), B.rows_, B.Data(), B.cols_, beta, C.Data(), B.cols_);
        }



};

template<>
struct BLAS<float>::func {
    constexpr static auto mul = cblas_sgemm;
    // constexpr static auto sum = cblas_saxpy;
};

template<>
struct BLAS<double>::func {
    constexpr static auto mul = cblas_dgemm;
    // constexpr static auto sum = cblas_daxpy;
};

} // namespace s21
