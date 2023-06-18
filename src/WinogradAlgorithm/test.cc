#include "Winograd.h"

#include "../sub/Utility/Utility.h"

using namespace s21;

typedef double fp_type;

int main() {
    int n = 32;
    Matrix<fp_type> A(n, n, [] { return Random::Normal<fp_type>(0, 15); });
    Matrix<fp_type> B(n, n, [] { return Random::Normal<fp_type>(0, 15); });
    // int k = 0;
    // Matrix<fp_type> A(n, n, [&] { return ++k; });
    // Matrix<fp_type> B(n, n, [&] { return ++k; });
    // Matrix<fp_type> A(n, n, [] { return Random::Int(1, 10); });
    // Matrix<fp_type> B(n, n, [] { return Random::Int(1, 10); });


    Matrix<fp_type> C2(n);


    // auto C1 = A * B;
    // Winograd<fp_type>::strassenWinogradMultiply(A, B, C2);
    // std::cout << C1 << '\n' << C2 << '\n';



    SStr::Print(Time::Compare<Time::ms>(10, [&] {
        auto C1 = A * B;
    }, [&] {
        Winograd<fp_type>::strassenWinogradMultiply(A, B, C2);
    }));

    return 0;
}
