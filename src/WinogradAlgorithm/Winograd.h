#pragma once

// #define MATRIX_1D

#include "../sub/matrix/matrix.h"

#include <vector>
#include <thread>

#include <memory>

// #include <mpi.h>

namespace s21a {

using namespace s21;

// Matrix(row, col)
// n m k
// C(n, m) = A(n, k) * B(k, m)


// Matrix(row, col)
// n m k
// C(n, m) = A(n, k) * B(k, m)

template<class T>
class Winograd {
    using M = Matrix<T>;
    using i_type = typename M::i_type;
    using data_t = typename M::base;

    public:

        // C(n, m) = A(n, k) * B(k, m)
        Winograd(i_type n);

        void Mul(const M &A, const M &B, M &C);

    private:
        struct Level {
            data_t A11, A12, A22, B11, B21, B22;
            data_t S1, S2, S3, S4, T1, T2, T3, T4;
            data_t R1, R2, R3, R4, R5, R6, R7;
            i_type n;

            Level(i_type n) : A11(n * n), A12(n * n), A22(n * n), B11(n * n), B21(n * n), B22(n * n),
                    S1(n * n), S2(n * n), S3(n * n), S4(n * n), T1(n * n), T2(n * n), T3(n * n), T4(n * n),
                    R1(n * n), R2(n * n), R3(n * n), R4(n * n), R5(n * n), R6(n * n), R7(n * n), n(n) {}
            
            virtual void MxABS(const T *A, const T *B) = 0;
            virtual void MxR(T *C) = 0;
        };

        struct LevelOdd : public Level {
            using Level::n;
            using Level::A11, Level::A12, Level::A22,
                Level::B11, Level::B21, Level::B22,
                Level::S1, Level::S2, Level::S3, Level::S4,
                Level::T1, Level::T2, Level::T3, Level::T4,
                Level::R1, Level::R2, Level::R3, Level::R4, Level::R5, Level::R6, Level::R7;

            LevelOdd(i_type n) : Level(n) {}
            void MxABS(const T *A, const T *B) override;
            void MxR(T *C) override;
        };

        struct LevelEven : public Level {
            using Level::n;
            using Level::A11, Level::A12, Level::A22,
                Level::B11, Level::B21, Level::B22,
                Level::S1, Level::S2, Level::S3, Level::S4,
                Level::T1, Level::T2, Level::T3, Level::T4,
                Level::R1, Level::R2, Level::R3, Level::R4, Level::R5, Level::R6, Level::R7;

            LevelEven(i_type n) : Level(n) {}
            void MxABS(const T *A, const T *B) override;
            void MxR(T *C) override;
        };

        std::vector<std::unique_ptr<Level>> L_;
        i_type n_, depth_;

        void sw_helper(const T *A, const T *B, T *C, int n, int l);
        void LastStep(const T *A, const T *B, T *C);

};

template<class T>
Winograd<T>::Winograd(i_type n) : n_(n) {
    // i_type f = std::cail(std::log2(n));

    if (n_ & (n_ - 1)) {
        // throw std::invalid_argument("n must be power of 2");
        bool even = (n % 2 == 0);
        if (!even) {
            ++n;
        }
        while ((n /= 2) != 1) {

            if (even)
                L_.push_back(std::make_unique<LevelEven>(n));
            else
                L_.push_back(std::make_unique<LevelOdd>(n));

            even = (n % 2 == 0);
            if (!even) {
                ++n;
            }
        }
    } else {
        while ((n /= 2) != 1) {
            L_.push_back(std::make_unique<LevelEven>(n));
        }
    }

    
}


template<class T>
void Winograd<T>::Mul(const M &A, const M &B, M &C) {
    if (A.GetCols() != n_ || B.GetCols() != n_ || C.GetCols() != n_ ||
        A.GetRows() != n_ || B.GetRows() != n_ || C.GetRows() != n_) {
        throw std::invalid_argument("Matrix size not match");
    }
    sw_helper(A.Data().data(), B.Data().data(), C.Data().data(), n_ + 1, 0);
}

template<class T>
void Winograd<T>::LevelEven::MxABS(const T *A, const T *B) {
    for (i_type i = 0, in = n; i < n; ++i, ++in) {
        for (i_type j = 0, jn = n; j < n; ++j, ++jn) {
            T a11 = A[i * n * 2 + j];
            T a12 = A[i * n * 2 + jn];
            T a21 = A[in * n * 2 + j];
            T a22 = A[in * n * 2 + jn];
            T b11 = B[i * n * 2 + j];
            T b12 = B[i * n * 2 + jn];
            T b21 = B[in * n * 2 + j];
            T b22 = B[in * n * 2 + jn];
            A11[i * n + j] = a11;
            A12[i * n + j] = a12;
            A22[i * n + j] = a22;
            B11[i * n + j] = b11;
            B21[i * n + j] = b21;
            B22[i * n + j] = b22;
            S1[i * n + j] = a21 + a22;
            S2[i * n + j] = a21 + a22 - a11;
            S3[i * n + j] = a11 - a21;
            S4[i * n + j] = a12 - a21 - a22 + a11;
            T1[i * n + j] = b12 - b11;
            T2[i * n + j] = b22 - b12 + b11;
            T3[i * n + j] = b22 - b12;
            T4[i * n + j] = b22 - b12 + b11 - b21;
        }
    }
}

template<class T>
void Winograd<T>::LevelEven::MxR(T *C) {
    for (i_type i = 0, in = n; i < n; ++i, ++in) {
        for (i_type j = 0, jn = n; j < n; ++j, ++jn) {
            T r1 = R1[i * n + j];
            T r7 = R7[i * n + j];
            T r16 = r1 + R6[i * n + j];
            T r165 = r16 + R5[i * n + j];
            C[i * n * 2 + j] = r1 + R2[i * n + j];
            C[i * n * 2 + jn] = r165 + R3[i * n + j];
            C[in * n * 2 + j] = r16 - R4[i * n + j] + r7;
            C[in * n * 2 + jn] = r165 + r7;
        }
    }
}

template<class T>
void Winograd<T>::LevelOdd::MxABS(const T *A, const T *B) {
    // std::cout << "LevelOdd::MxABS n = " << n << "\n";
    for (i_type i = 0, in = n; i < n; ++i, ++in) {
        for (i_type j = 0, jn = n; j < n; ++j, ++jn) {
            i_type cn = n * 2 - 1;
            T a11 = A[i * cn + j];
            T a12 = (jn == n * 2 - 1) ? 0 : A[i * cn + jn];
            T a21 = (in == n * 2 - 1) ? 0 : A[in * cn + j];
            T a22 = (in == n * 2 - 1 || jn == n * 2 - 1) ? 0 : A[in * cn + jn];
            T b11 = B[i * cn + j];
            T b12 = (jn == n * 2 - 1) ? 0 : B[i * cn + jn];
            T b21 = (in == n * 2 - 1) ? 0 : B[in * cn + j];
            T b22 = (in == n * 2 - 1 || jn == n * 2 - 1) ? 0 : B[in * cn + jn];
            A11[i * n + j] = a11;
            A12[i * n + j] = a12;
            A22[i * n + j] = a22;
            B11[i * n + j] = b11;
            B21[i * n + j] = b21;
            B22[i * n + j] = b22;
            S1[i * n + j] = a21 + a22;
            S2[i * n + j] = a21 + a22 - a11;
            S3[i * n + j] = a11 - a21;
            S4[i * n + j] = a12 - a21 - a22 + a11;
            T1[i * n + j] = b12 - b11;
            T2[i * n + j] = b22 - b12 + b11;
            T3[i * n + j] = b22 - b12;
            T4[i * n + j] = b22 - b12 + b11 - b21;
        }
    }
}

template<class T>
void Winograd<T>::LevelOdd::MxR(T *C) {
    Matrix<T> MxR(n * 2);
    for (i_type i = 0, in = n; i < n; ++i, ++in) {
        for (i_type j = 0, jn = n; j < n; ++j, ++jn) {
            T r1 = R1[i * n + j];
            T r7 = R7[i * n + j];
            T r16 = r1 + R6[i * n + j];
            T r165 = r16 + R5[i * n + j];

            i_type cn = n * 2 - 1;
            C[i * cn + j] = r1 + R2[i * n + j];
            if (!(jn == n * 2 - 1))
                C[i * cn + jn] = r165 + R3[i * n + j];
            if (!(in == n * 2 - 1))
                C[in * cn + j] = r16 - R4[i * n + j] + r7;
            if (!(in == n * 2 - 1) && !(jn == n * 2 - 1))
                C[in * cn + jn] = r165 + r7;
        }
    }
}

template<class T>
void Winograd<T>::LastStep(const T *A, const T *B, T *C) {
    T a11 = A[0];
    T a12 = A[1];
    T a21 = A[2];
    T a22 = A[3];
    T b11 = B[0];
    T b12 = B[1];
    T b21 = B[2];
    T b22 = B[3];
    C[0] = a11 * b11 + a12 * b21;
    C[1] = a11 * b12 + a12 * b22;
    C[2] = a21 * b11 + a22 * b21;
    C[3] = a21 * b12 + a22 * b22;
}

template<class T>
void Winograd<T>::sw_helper(const T *A, const T *B, T *C, int n, int l) {

    n /= 2;

    if (n == 1) {
        T a11 = A[0];
        T a12 = A[1];
        T a21 = A[2];
        T a22 = A[3];
        T b11 = B[0];
        T b12 = B[1];
        T b21 = B[2];
        T b22 = B[3];
        C[0] = a11 * b11 + a12 * b21;
        C[1] = a11 * b12 + a12 * b22;
        C[2] = a21 * b11 + a22 * b21;
        C[3] = a21 * b12 + a22 * b22;
        return;
    }

    auto &L = *L_[l++].get();

    L.MxABS(A, B);

    sw_helper(L.A11.data(), L.B11.data(), L.R1.data(), n + 1, l);
    sw_helper(L.A12.data(), L.B21.data(), L.R2.data(), n + 1, l);
    sw_helper(L.S4.data(), L.B22.data(), L.R3.data(), n + 1, l);
    sw_helper(L.A22.data(), L.T4.data(), L.R4.data(), n + 1, l);
    sw_helper(L.S1.data(), L.T1.data(), L.R5.data(), n + 1, l);
    sw_helper(L.S2.data(), L.T2.data(), L.R6.data(), n + 1, l);
    sw_helper(L.S3.data(), L.T3.data(), L.R7.data(), n + 1, l);

    L.MxR(C);
}


} // namespace s21

