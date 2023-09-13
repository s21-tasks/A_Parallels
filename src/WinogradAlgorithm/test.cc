

// #include "Winograd.h"
// #include "staticMatrix.h"
#include "Winograd.h"

#include "Winograd copy 2.h"

#include "Winograd copy_func_ar.h"


#include "../sub/utility/utility.h"

using namespace s21;

typedef int fp_type;

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
    // Matrix<fp_type> C5(n, n);

    s21a::Winograd<fp_type> W2(n);
    Winograd<fp_type> W(n);
    s21b::Winograd<fp_type> W3(n);


    std::cout << "START time test\n";
    auto T = Time::Now();
    W3.Mul(A, B, C2);
    std::cout << Time::Duration<Time::ms>(T) << " ms\n";
    T = Time::Now();
    auto C1 = A * B;
    std::cout << Time::Duration<Time::ms>(T) << " ms\n";
    T = Time::Now();
    W2.Mul(A, B, C3);
    std::cout << Time::Duration<Time::ms>(T) << " ms\n";
    

    // std::cout << C1 << '\n' << C2 << '\n';
    std::cout << "(C1 == C2) = " << (C1 == C2) << '\n';
    std::cout << "(C1 == C3) = " << (C1 == C3) << '\n';
    // std::cout << (C1 == C3) << '\n';

    // SStr::Print(Time::Compare<Time::ms>(1, [&] {
    //     auto C1 = A * B;
    // }, [&] {
    //     W.MulAdj();
    // }, [&] {
    //     // Winograd<fp_type>::MulOneTime(A, B, C2);
    // }, [&] {
    //     // Winograd<fp_type>::MulThreadM(A, B, C4);
    // }));

    return 0;
}
