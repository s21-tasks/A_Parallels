

// #include "Winograd.h"
// #include "staticMatrix.h"
#include "Winograd.h"

#include "../sub/utility/utility.h"

#include "../sub/matrix/blas.h"

using namespace s21;

typedef float fp_type;

int main() {

    // std::cout << std::ceil(std::log2(257)) << '\n';
    // int l = 1, m = 2, k = 3;
    // std::cout << std::max({l, m, k}) << '\n';

    // int k = 16;
    // while ((k /= 2) != 1) {
    //     std::cout << k << '\n';
    // }

    int n = 1024;
    Matrix<fp_type> A(n, n, [] { return Random::Int(-10, 10); });
    Matrix<fp_type> B(n, n, [] { return Random::Int(-10, 10); });
    // Matrix<fp_type> A(n, n, 5);
    // Matrix<fp_type> B(n, n, 5);
    // A.Resize(8, 8, 0);
    // B.Resize(8, 8, 0);
    // n = 8;


    Matrix<fp_type> C2(n, n);
    Matrix<fp_type> C3(n, n);
    Matrix<fp_type> C4(n, n);
    Matrix<fp_type> C5(n, n);

    Winograd<fp_type> W(n);


    std::cout << "START time test\n";
    auto T = Time::Now();
    W.Mul(A, B, C2);
    std::cout << Time::Duration<Time::ms>(T) << " ms\n";
    T = Time::Now();
    // auto C1 = A * B;
    Blas<fp_type>::Mul(A, B, C5);
    std::cout << Time::Duration<Time::ms>(T) << " ms\n";
    T = Time::Now();
    // W2.Mul(A, B, C3);
    Matrix<fp_type>::Mul(A, B, C3);
    std::cout << Time::Duration<Time::ms>(T) << " ms\n";
    

    // std::cout << C5 << '\n' << C2 << '\n';
    std::cout << "(C5 == C2) = " << (C5 == C2) << '\n';
    std::cout << "(C5 == C3) = " << (C5 == C3) << '\n';
    std::cout << "(C2 == C3) = " << (C2 == C3) << '\n';




    // Matrix<fp_type> C9(n, n);
    // Matrix<fp_type> C8(n, n);

    // Blas<fp_type>::MulAdj(A, B, C9);
    // Matrix<fp_type>::MulAdj(A, B, C8);
    // std::cout << C9 << ' ' << C8 << '\n';
    // std::cout << (C9 == C8) << '\n';

    return 0;
}
