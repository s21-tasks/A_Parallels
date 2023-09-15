#define LEVEL_LOAD_TEST__P

#include "../winograd.h"
#include "../../sub/utility/utility.h"

using namespace s21;

template<class T>
T RS() {
    return Random::Easy<T>::F(static_cast<T>(-10), static_cast<T>(10));
}

template<class T>
void F(unsigned int N) {
    static int counter;
    using W = typename Basic::Winograd<T>;
    W::Level22::load = 0;
    W::LevelClassic::load = 0;
    W::LevelEven::level_load.clear();
    W::LevelOdd::level_load.clear();
    Matrix<T> A(N, N, RS<T>);
    Matrix<T> B(N, N, RS<T>);
    Matrix<T> C(N, N);
    for (int k = 0; k < 10; ++ k) {
        W::Mul(A, B, C);
    }
    std::cout << "\033[1;32m";
    std::cout << "Test Level Load " << ++counter << " {N = " << N << ", Unit = ms}:\n";
    std::cout << "\033[1;37m";
    std::cout << "\tLevel Even: {";
    for (auto &p : W::LevelEven::level_load) {
        std::cout << "\"" << p.first << " - " << p.second / 10000000 << "\" ";
    }
    std::cout << "}\n\tLevel Odd: {";
    for (auto &p : W::LevelOdd::level_load) {
        std::cout << "\"" << p.first << " - " << p.second / 10000000 << "\" ";
    }
    std::cout << "}\n\tLevel 2x2: " << W::Level22::load / 10000000;
    std::cout << "\n\tLevel Classic: " << W::LevelClassic::load / 10000000 << "\n\n";
}

int main() {
    F<float>(1024);
    F<float>(2048);
    F<float>(1000);
    F<float>(1011);
    F<float>(1012);
    F<float>(1013);
    F<float>(999);

    return 0;
}
