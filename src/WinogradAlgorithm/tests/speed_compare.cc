#include "tests.h"
#include "../../sub/matrix/blas.h"

#include <fstream>
#include <map>

using namespace s21;

template<class T, class Unit>
void All(unsigned int N, unsigned int repeat) {
    static int counter = 0;

    Matrix<T> A(N, N, RS<T>);
    Matrix<T> B(N, N, RS<T>);
    Matrix<T> C(N, N);
    auto result = Time::Compare<Unit>(repeat, [&] {
            Basic::Winograd<T>::Mul(A, B, C);
        }, [&] {
            Parallel::Winograd<T>::Mul(A, B, C);
        }, [&] {
            BLAS<T>::Mul(A, B, C);
        }, [&] {
            if (N <= 512) Matrix<T>::Mul(A, B, C);
        }, [&] {
            Matrix<T>::MulParallel(A, B, C, 4);
        }
    );
    std::cout << "\033[1;32m";
    std::cout << "\nTest all algorithms " << ++counter << " {N = " << N << ". Type = " << typeid(T).name() << "}:\n";
    std::cout << "\033[1;37m";
    std::cout << '\t' << result[0] << ' ' << Time::Prefix<Unit>::val << "\t - Basic Winograd\n";
    std::cout << '\t' << result[1] << ' ' << Time::Prefix<Unit>::val << "\t - Parallel Winograd\n";
    std::cout << '\t' << result[2] << ' ' << Time::Prefix<Unit>::val << "\t - Blas\n";
    if (N <= 512)
        std::cout << '\t' << result[3] << ' ' << Time::Prefix<Unit>::val << "\t - Classic\n";
    else
        std::cout << '\t' << "дохуя\t - Classic\n";
    std::cout << '\t' << result[4] << ' ' << Time::Prefix<Unit>::val << "\t - Classic Parallel 4\n";
}



// odd cap and winograd cap try to find the best value
template<class T, class Unit>
void Cap(unsigned int N, unsigned int repeat) {

    static int counter = 0;
    Matrix<T> A(N, N, RS<T>);
    Matrix<T> B(N, N, RS<T>);
    Matrix<T> C(N, N);
    
    std::vector<unsigned int> caps{15, 16, 25, 33, 34, 47, 63, 64, 81, 127, 128, 129};

    for (unsigned int k = 0; k < caps.size(); ++k) {
        for (unsigned int g = k; g < caps.size(); ++g) {
            auto result = Time::Compare<Unit>(repeat, [&] {
                    Basic::Winograd<T> w(N, caps[k], caps[g]);
                    w.Mul(A, B, C, caps[g], caps[k]);
                }, [&] {
                    Parallel::Winograd<T> w(N, caps[k], caps[g]);
                    w.Mul(A, B, C, caps[g], caps[k]);
                }, [&] {
                    BLAS<T>::Mul(A, B, C);
                }
            );
            std::cout << "\033[1;32m";
            std::cout << "\nTest cap " << ++counter << " {N = " << N << ". Type = " << typeid(T).name() << "}:\n";
            std::cout << "\033[1;37m";
            std::cout << "\t\tOdd cap = " << caps[g] << "\tWinograd cap = " << caps[k] << '\n';
            std::cout << '\t' << result[0] << ' ' << Time::Prefix<Unit>::name << "\t - Basic Winograd\n";
            std::cout << '\t' << result[1] << ' ' << Time::Prefix<Unit>::name << "\t - Parallel Winograd\n";
            std::cout << '\t' << result[2] << ' ' << Time::Prefix<Unit>::name << "\t - Blas\n";
        }
    }
}

// template<class T>
// void ToFile3(std::ofstream &file, const typename Matrix<T>::i_type N, int repeat) {
//     Matrix<T> A(N, N, Random::Easy<T>::F(static_cast<T>(-10), static_cast<T>(10)));
//     Matrix<T> B(N, N, Random::Easy<T>::F(static_cast<T>(-10), static_cast<T>(10)));
//     Matrix<T> C1(N, N);
//     Matrix<T> C2(N, N);
//     Matrix<T> C3(N, N);
//     auto result = Time::Compare<Time::ns>(repeat, [&] {
//         Winograd<T>::Mul(A, B, C1);
//     }, [&] {
//         Blas<T>::Mul(A, B, C2);
//     }, [&] {
//         Matrix<T>::Mul(A, B, C3);
//     });
//     file << N << ' ' << result[0] << ' ' << result[1] << ' ' << result[2] << '\n';
// }

// template<class T>
// void ToFile2(std::ofstream &file, const typename Matrix<T>::i_type N, int repeat) {
//     Matrix<T> A(N, N, Random::Easy<T>::F(static_cast<T>(-10), static_cast<T>(10)));
//     Matrix<T> B(N, N, Random::Easy<T>::F(static_cast<T>(-10), static_cast<T>(10)));
//     Matrix<T> C1(N, N);
//     Matrix<T> C2(N, N);
//     auto result = Time::Compare<Time::ns>(repeat, [&] {
//         Winograd<T>::Mul(A, B, C1);
//     }, [&] {
//         Blas<T>::Mul(A, B, C2);
//     });
//     file << N << ' ' << result[0] << ' ' << result[1] << '\n';
// }

int main() {

    All<float, Time::ns>(8, 100000);
    All<float, Time::ns>(25, 10000);
    All<float, Time::mcs>(64, 10000);
    // All<double, Time::mcs>(64, 10000);
    All<float, Time::mcs>(78, 1000);
    All<float, Time::mcs>(128, 1000);
    All<float, Time::mcs>(133, 1000);
    All<float, Time::mcs>(211, 1000);
    All<float, Time::mcs>(255, 400);
    All<float, Time::mcs>(256, 300);
    All<float, Time::ms>(512, 100);
    All<float, Time::ms>(1024, 6);
    // All<double, Time::ms>(1024, 6);
    // All<float, Time::ms>(2048, 3);

    // // Cap<float, Time::ns>(8, 100000);
    // // Cap<float, Time::ns>(25, 10000);
    // // Cap<float, Time::mcs>(64, 10000);
    // // Cap<double, Time::mcs>(64, 10000);
    // Cap<float, Time::mcs>(78, 1000);
    // Cap<float, Time::mcs>(128, 1000);
    // // Cap<float, Time::mcs>(133, 1000);
    // Cap<float, Time::mcs>(211, 1000);
    // Cap<float, Time::mcs>(255, 400);
    // Cap<float, Time::mcs>(256, 300);
    // // Cap<float, Time::ms>(512, 100);
    // // Cap<float, Time::ms>(1024, 6);
    // // Cap<double, Time::ms>(1024, 6);
    // // Cap<float, Time::ms>(2048, 3);

    // G<float, Time::sec>(4096, 2);


    // std::ofstream file("speed_test_data_2_512.txt");
    // if (!file.is_open()) {
    //     std::cout << "File not found\n";
    //     return 0;
    // }
    // file.clear();
    // auto time_point = Time::Now();
    // for (int n = 2, r = 10000; n < 128; ++n, r -= 50) {
    //     ToFile3<double>(file, n, r);
    // }
    // for (int n = 128, r = 2000; n < 512; ++n, r -= 4) {
    //     ToFile3<double>(file, n, 1000);
    // }
    // // for (int n = 512, r = 200; n < 1024; n += 7, r -= 2) {
    // //     ToFile3<double>(file, n, 1000);
    // // }
    // // ToFile3<double>(file, 1024, 8);
    // // int c = 0;
    // // for (int n = 1025; n < 2048; ) {
    // //     ToFile2<double>(file, n, 5);
    // //     if (++c % 2 == 0) {
    // //         n += 7;
    // //     } else {
    // //         n += 8;
    // //     }
    // // }
    // // ToFile3<double>(file, 2048, 5);
    // // for (int n = 2048; n < 4096; n += 31) {
    // //     ToFile2<double>(file, n, 3);
    // // }
    // // ToFile2<double>(file, 4096, 5);
    // file.close();
    // std::cout << "Time: " << Time::Duration<Time::sec>(time_point) << " sec\n";

    return 0;
}
