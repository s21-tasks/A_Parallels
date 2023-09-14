#pragma once

// #define MATRIX_1D

#include "../sub/matrix/matrix.h"

#include <vector>
#include <thread>
#include <memory>

// #include <mpi.h>

namespace s21 {

using namespace s21;

template<class T>
class Winograd {
    using M = Matrix<T>;
    using i_type = typename M::i_type;
    using data_t = typename M::base;

    public:
        Winograd(i_type n, i_type square_cap = 128);
        void Mul(const M &A, const M &B, M &C);
        static void Mul(const M &A, const M &B, M &C);

    private:
        struct Level {
            virtual void SW(const T *A, const T *B, T *C) = 0;
            Level *next;
            virtual ~Level() = default;
        };

        struct LevelEven : public Level {
            using Level::next;

            T *A11, *A12, *A22, *B11, *B21, *B22;
            T *S1, *S2, *S3, *S4, *T1, *T2, *T3, *T4;
            T *R1, *R2, *R3, *R4, *R5, *R6, *R7;
            i_type n;

            LevelEven(i_type n) : n(n) {
                A11 = new T[n * n];
                A12 = new T[n * n];
                A22 = new T[n * n];
                B11 = new T[n * n];
                B21 = new T[n * n];
                B22 = new T[n * n];
                S1 = new T[n * n];
                S2 = new T[n * n];
                S3 = new T[n * n];
                S4 = new T[n * n];
                T1 = new T[n * n];
                T2 = new T[n * n];
                T3 = new T[n * n];
                T4 = new T[n * n];
                R1 = new T[n * n];
                R2 = new T[n * n];
                R3 = new T[n * n];
                R4 = new T[n * n];
                R5 = new T[n * n];
                R6 = new T[n * n];
                R7 = new T[n * n];
            }

            virtual ~LevelEven() {
                delete[] A11;
                delete[] A12;
                delete[] A22;
                delete[] B11;
                delete[] B21;
                delete[] B22;
                delete[] S1;
                delete[] S2;
                delete[] S3;
                delete[] S4;
                delete[] T1;
                delete[] T2;
                delete[] T3;
                delete[] T4;
                delete[] R1;
                delete[] R2;
                delete[] R3;
                delete[] R4;
                delete[] R5;
                delete[] R6;
                delete[] R7;
            }

            virtual void SW(const T *A, const T *B, T *C) override;
        };

        struct LevelOdd final : public LevelEven {
            using LP = LevelEven;
            using LP::A11, LP::A12, LP::A22,
                LP::B11, LP::B21, LP::B22,
                LP::S1, LP::S2, LP::S3, LP::S4,
                LP::T1, LP::T2, LP::T3, LP::T4,
                LP::R1, LP::R2, LP::R3, LP::R4, LP::R5, LP::R6, LP::R7;
            using LP::n, LP::next;

            LevelOdd(i_type n) : LevelEven(n) {}

            void SW(const T *A, const T *B, T *C) override;
        };

        struct LevelClassic final : public Level {
            i_type n;
            using Level::next;

            LevelClassic(i_type n) : n(n) {}

            void SW(const T *A, const T *B, T *C) override;
        };

        struct Level22 final : public Level {
            using Level::next;

            void SW(const T *A, const T *B, T *C) override;   
        };
        
        std::vector<std::unique_ptr<Level>> L_;

        i_type n_;

};

template
struct Level {
    virtual void SW(const T *A, const T *B, T *C) = 0;
    Level *next;
    virtual ~Level() = default;
};

template<class T>
Winograd<T>::Winograd(i_type n, i_type square_cap) : n_(n) {
    bool even = (n % 2 == 0);
    if (!even) {
        ++n;
    }
    bool cap = false;

    while ((n /= 2) != 1) {
        if (even)
            L_.push_back(std::make_unique<LevelEven>(n));
        else
            L_.push_back(std::make_unique<LevelOdd>(n));

        even = (n % 2 == 0);
        if (!even) {
            if (n < square_cap) {
                std::cout << n << " cap\n";
                L_.push_back(std::make_unique<LevelClassic>(n));
                cap = true;
                break;
            }
            ++n;
        }
    }

    if (!cap) {
        L_.push_back(std::make_unique<Level22>());
    }

    for (int i = L_.size() - 1; i > 0; --i) {
        L_[i - 1]->next = &(*L_[i]);
    }

}


template<class T>
void Winograd<T>::Mul(const M &A, const M &B, M &C) {
    if (A.GetCols() != n_ || B.GetCols() != n_ || C.GetCols() != n_ ||
        A.GetRows() != n_ || B.GetRows() != n_ || C.GetRows() != n_) {
        throw std::invalid_argument("Matrix size not match");
    }
    L_[0]->SW(A.Data().data(), B.Data().data(), C.Data().data());
}

template<class T>
void Winograd<T>::Level22::SW(const T *A, const T *B, T *C) {
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
void Winograd<T>::LevelClassic::SW(const T *A, const T *B, T *C) {
    for (i_type i = 0; i < n; ++i) {
        for (i_type j = 0; j < n; ++j) {
            T sum = 0;
            for (i_type k = 0; k < n; ++k) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
};

template<class T>
void Winograd<T>::LevelEven::SW(const T *A, const T *B, T *C) {

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


    next->SW(A11, B11, R1);
    next->SW(A12, B21, R2);
    next->SW(S4, B22, R3);
    next->SW(A22, T4, R4);
    next->SW(S1, T1, R5);
    next->SW(S2, T2, R6);
    next->SW(S3, T3, R7);


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
void Winograd<T>::LevelOdd::SW(const T *A, const T *B, T *C) {

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


    next->SW(A11, B11, R1);
    next->SW(A12, B21, R2);
    next->SW(S4, B22, R3);
    next->SW(A22, T4, R4);
    next->SW(S1, T1, R5);
    next->SW(S2, T2, R6);
    next->SW(S3, T3, R7);


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


} // namespace s21

