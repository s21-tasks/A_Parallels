#pragma once

#include "../sub/Matrix/Matrix.h"

#include <vector>

namespace s21 {

template<class T>
struct WMatrix {
    const Matrix<T> &m_;
    int c1_, c2_, r1_, r2_;
    WMatrix(const Matrix<T> &m, int c1, int c2, int r1, int r2) :
        m_(m), c1_(c1), c2_(c2), r1_(r1), r2_(r2) {}
    WMatrix(const Matrix<T> &m) :
        m_(m), c1_(0), c2_(m.GetCols()), r1_(0), r2_(m.GetRows()) {}

    const T &operator()(int r, int c) const {
        return m_(r + r1_, c + c1_);
    }
};

template<class T>
class Winograd {

    public:
        static Matrix<T> Mul(Matrix<T> &A, Matrix<T> &B) {
            int m = A.GetRows();
            int n = A.GetCols();
            int k = B.GetCols();

            int d = n / 2;


            Matrix<T> C(m, k, 0);

            std::vector<T> row_f(m, 0.0);
            std::vector<T> col_f(k, 0.0);

            for (int i = 0; i < m; ++i) {
                // row_f[i] = 0;
                for (int j = 0; j < d; ++j) {
                    row_f[i] += A(i, 2 * j) * A(i, 2 * j + 1);
                }
            }

            for (int i = 0; i < k; ++i) {
                // clo_f[i] = 0;
                for (int j = 0; j < d; ++j) {
                    col_f[i] += B(2 * j, i) * B(2 * j + 1, i);
                }
            }

            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < k; ++j) {
                    C(i, j) = -1.0 * (row_f[i] + col_f[j]);
                    for (int k = 0; k < d; ++k) {
                        C(i, j) += (A(i, 2 * k) + B(2 * k + 1, j)) * (A(i, 2 * k + 1) + B(2 * k, j));
                        if (n % 2 !=0) {
                            C(i, j) += A(i, n - 1) * B(n - 1, j);
                        }
                    }
                }
            }

        

            return C;

        }

        static void strassenWinogradMultiply(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C);

        static void strassenWinogradMultiply(const WMatrix<T> &A, const WMatrix<T> &B, Matrix<T> &C, int &count);
};



template<class T>
void Winograd<T>::strassenWinogradMultiply(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C) {

    // WMatrix<T> WA(A, 0, A.GetCols() - 1, 0, A.GetRows() - 1);
    // WMatrix<T> BA(B, 0, B.GetCols() - 1, 0, B.GetRows() - 1);
    int count = 0;
    strassenWinogradMultiply(WMatrix<T>(A), WMatrix<T>(B), C, count);
    std::cout << "count: " << count << '\n';

}

template<class T>
void Winograd<T>::strassenWinogradMultiply(const WMatrix<T> &A, const WMatrix<T> &B, Matrix<T> &C, int &count) {
    int n = A.r2_ - A.r1_;

    // std::cout << A.c1_ << ' ' << A.c2_ << ' ' << A.r1_ << ' ' << A.r2_ << ' ' << B.c1_ << ' ' << B.c2_ << ' ' << B.r1_ << ' ' << B.r2_ << '\n';

    if (n == 1) {
        C(0, 0) = A(0, 0) * B(0, 0);
        ++count;
        return;
    }
    // if (n < 1) {
    //     std::cout << A.c1_ << ' ' << A.c2_ << ' ' << A.r1_ << ' ' << A.r2_ << ' ' << B.c1_ << ' ' << B.c2_ << ' ' << B.r1_ << ' ' << B.r2_ << '\n';
    //     return;
    // }

    int newSize = n / 2;

    Matrix<T> S1(newSize);
    Matrix<T> S2(newSize);
    Matrix<T> S3(newSize);
    Matrix<T> S4(newSize);

    Matrix<T> T1(newSize);
    Matrix<T> T2(newSize);
    Matrix<T> T3(newSize);
    Matrix<T> T4(newSize);

    int newAc = (A.c2_ + A.c1_) / 2;
    int newAr = (A.r2_ + A.r1_) / 2;
    int newBc = (B.c2_ + B.c1_) / 2;
    int newBr = (B.r2_ + B.r1_) / 2;

    WMatrix<T> A11(A.m_, A.c1_, newAc, A.r1_, newAr);
    WMatrix<T> A12(A.m_, newAc, A.c2_, A.r1_, newAr);
    WMatrix<T> A21(A.m_, A.c1_, newAc, newAr, A.r2_);
    WMatrix<T> A22(A.m_, newAc, A.c2_, newAr, A.r2_);

    WMatrix<T> B11(B.m_, B.c1_, newBc, B.r1_, newBr);
    WMatrix<T> B12(B.m_, newBc, B.c2_, B.r1_, newBr);
    WMatrix<T> B21(B.m_, B.c1_, newBc, newBr, B.r2_);
    WMatrix<T> B22(B.m_, newBc, B.c2_, newBr, B.r2_);

    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {

            // int ai = i + A.r1_, aj = j + A.c1_, bi = B.r1_, bj = B.c1_;

            S1(i, j) = A21(i, j) + A22(i, j);
            S2(i, j) = S1(i, j) - A11(i, j);
            S3(i, j) = A11(i, j) - A21(i, j);
            S4(i, j) = A12(i, j) - S2(i, j);

            T1(i, j) = B12(i, j) - B11(i, j);
            T2(i, j) = B22(i, j) - T1(i, j);
            T3(i, j) = B22(i, j) - B12(i, j);
            T4(i, j) = T2(i, j) - B21(i, j);

            count += 8;

        }
    }

    Matrix<T> R1(newSize);
    Matrix<T> R2(newSize);
    Matrix<T> R3(newSize);
    Matrix<T> R4(newSize);
    Matrix<T> R5(newSize);
    Matrix<T> R6(newSize);
    Matrix<T> R7(newSize);

    strassenWinogradMultiply(A11, B11, R1, count);
    strassenWinogradMultiply(A12, B21, R2, count);
    strassenWinogradMultiply(WMatrix<T>(S4), B22, R3, count);
    strassenWinogradMultiply(A22, WMatrix<T>(T4), R4, count);
    strassenWinogradMultiply(WMatrix<T>(S1), WMatrix<T>(T1), R5, count);
    strassenWinogradMultiply(WMatrix<T>(S2), WMatrix<T>(T2), R6, count);
    strassenWinogradMultiply(WMatrix<T>(S3), WMatrix<T>(T3), R7, count);

    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            T r16 = R1(i, j) + R6(i, j);
            T r165 = r16 + R5(i, j);

            C(i, j) = R1(i, j) + R2(i, j);
            C(i, j + newSize) = r165 + R3(i, j);
            C(i + newSize, j) = r16 - R4(i, j) + R7(i, j);
            C(i + newSize, j + newSize) = r165 + R7(i, j);

            count += 4;
        }
    }
}







// template<class T>
// void Winograd<T>::strassenWinogradMultiply(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C) {
//     int n = A.GetRows();
//     if (n == 1) {
//         C[0][0] = A[0][0] * B[0][0];
//         return;
//     }
//     int newSize = n / 2;
//     Matrix<T> A11(newSize);
//     Matrix<T> A12(newSize);
//     // Matrix<T> A21(newSize);
//     Matrix<T> A22(newSize);
//     Matrix<T> B11(newSize);
//     // Matrix<T> B12(newSize);
//     Matrix<T> B21(newSize);
//     Matrix<T> B22(newSize);
//     Matrix<T> S1(newSize);
//     Matrix<T> S2(newSize);
//     Matrix<T> S3(newSize);
//     Matrix<T> S4(newSize);
//     Matrix<T> T1(newSize);
//     Matrix<T> T2(newSize);
//     Matrix<T> T3(newSize);
//     Matrix<T> T4(newSize);
//     for (int i = 0; i < newSize; i++) {
//         for (int j = 0; j < newSize; j++) {
//             A11(i, j) = A(i, j);
//             A12(i, j) = A(i, j + newSize);
//             // A21(i, j) = A(i + newSize, j); //
//             A22(i, j) = A(i + newSize, j + newSize);
//             B11(i, j) = B(i, j);
//             // B12(i, j) = B(i, j + newSize); //
//             B21(i, j) = B(i + newSize, j);
//             B22(i, j) = B(i + newSize, j + newSize);
//             S1(i, j) = A(i + newSize, j) + A22(i, j);
//             S2(i, j) = S1(i, j) - A11(i, j);
//             S3(i, j) = A11(i, j) - A(i + newSize, j);
//             S4(i, j) = A12(i, j) - S2(i, j);
//             T1(i, j) = B(i, j + newSize) - B11(i, j);
//             T2(i, j) = B22(i, j) - T1(i, j);
//             T3(i, j) = B22(i, j) - B(i, j + newSize);
//             T4(i, j) = T2(i, j) - B21(i, j);
//         }
//     }
//     Matrix<T> R1(newSize);
//     Matrix<T> R2(newSize);
//     Matrix<T> R3(newSize);
//     Matrix<T> R4(newSize);
//     Matrix<T> R5(newSize);
//     Matrix<T> R6(newSize);
//     Matrix<T> R7(newSize);
//     strassenWinogradMultiply(A11, B11, R1);
//     strassenWinogradMultiply(A12, B21, R2);
//     strassenWinogradMultiply(S4, B22, R3);
//     strassenWinogradMultiply(A22, T4, R4);
//     strassenWinogradMultiply(S1, T1, R5);
//     strassenWinogradMultiply(S2, T2, R6);
//     strassenWinogradMultiply(S3, T3, R7);
//     for (int i = 0; i < newSize; i++) {
//         for (int j = 0; j < newSize; j++) {
//             T r16 = R1(i, j) + R6(i, j);
//             T r165 = r16 + R5(i, j);
//             C(i, j) = R1(i, j) + R2(i, j);
//             C(i, j + newSize) = r165 + R3(i, j);
//             C(i + newSize, j) = r16 - R4(i, j) + R7(i, j);
//             C(i + newSize, j + newSize) = r165 + R7(i, j);
//         }
//     }
// }

} // namespace s21


