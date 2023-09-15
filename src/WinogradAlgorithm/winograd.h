#pragma once

#include "../sub/matrix/matrix.h"

#include <vector>
#include <thread>
#include <memory>

#ifdef LEVEL_LOAD_TEST__P
    #include <map>
    #include "../sub/utility/m_time.h"
#endif

namespace s21 {

using namespace s21;

namespace Basic {

template<class T>
class Winograd {
    using M = Matrix<T>;
    using i_type = typename M::i_type;
    using data_t = typename M::base;

    public:
        Winograd(i_type n, i_type odd_cap = 25, i_type winograd_cap = 15);
        void Execute(const M &A, const M &B, M &C);
        static void Mul(const M &A, const M &B, M &C,
                        i_type odd_cap = 25, i_type winograd_cap = 15);
        
        virtual ~Winograd() = default;

#ifdef LEVEL_LOAD_TEST__P
    public:
#else
    protected:
#endif
        Winograd() = default;

        struct Level;
        struct LevelEven;
        struct LevelOdd;

        struct LevelClassic;
        struct Level22;

        struct PointerBase;
        
        std::vector<std::unique_ptr<Level>> L_;

        i_type n_;
};

template<class T>
void Winograd<T>::Mul(const M &A, const M &B, M &C,
            i_type odd_cap, i_type winograd_cap) {

    i_type n = A.GetCols();
    if (B.GetCols() != n || C.GetCols() != n ||
        A.GetRows() != n || B.GetRows() != n || C.GetRows() != n) {
        throw std::invalid_argument("Matrix size not match");
    }
    Winograd<T> W(n, odd_cap, winograd_cap);
    W.Execute(A, B, C);
}

template<class T>
struct Winograd<T>::Level {
    virtual void SW(const T *A, const T *B, T *C) = 0;
    Level *next;
    virtual ~Level() = default;
};

template<class T>
struct Winograd<T>::PointerBase {
    T *A11, *A12, *A22, *B11, *B21, *B22;
    T *S1, *S2, *S3, *S4, *T1, *T2, *T3, *T4;
    T *R1, *R2, *R3, *R4, *R5, *R6, *R7;
    i_type n;

    PointerBase(i_type n);
    virtual ~PointerBase();
};

template<class T>
struct Winograd<T>::LevelEven : public Level, public PointerBase {
    using Level::next;
    using PB = PointerBase;
    using PB::n;
    using PB::A11, PB::A12, PB::A22, PB::B11, PB::B21, PB::B22,
        PB::S1, PB::S2, PB::S3, PB::S4, PB::T1, PB::T2, PB::T3, PB::T4,
        PB::R1, PB::R2, PB::R3, PB::R4, PB::R5, PB::R6, PB::R7;

#ifdef LEVEL_LOAD_TEST__P
    inline static std::map<i_type, int64_t> level_load;
#endif

    LevelEven(i_type n) : PointerBase(n) {}
    virtual ~LevelEven() = default;
    virtual void SW(const T *A, const T *B, T *C) override;
};

template<class T>
struct Winograd<T>::LevelOdd : public Level, public PointerBase {
    using Level::next;
    using PB = PointerBase;
    using PB::n;
    using PB::A11, PB::A12, PB::A22, PB::B11, PB::B21, PB::B22,
        PB::S1, PB::S2, PB::S3, PB::S4, PB::T1, PB::T2, PB::T3, PB::T4,
        PB::R1, PB::R2, PB::R3, PB::R4, PB::R5, PB::R6, PB::R7;
    
#ifdef LEVEL_LOAD_TEST__P
    inline static std::map<i_type, int64_t> level_load;
#endif

    LevelOdd(i_type n) : PointerBase(n) {}
    virtual ~LevelOdd() = default;
    virtual void SW(const T *A, const T *B, T *C) override;
};

template<class T>
struct Winograd<T>::LevelClassic final : public Level {
    i_type n;

#ifdef LEVEL_LOAD_TEST__P
    inline static int64_t load;
#endif

    LevelClassic(i_type n) : n(n) {}
    void SW(const T *A, const T *B, T *C) override;
};

template<class T>
struct Winograd<T>::Level22 final : public Level {
    using Level::next;

#ifdef LEVEL_LOAD_TEST__P
    inline static int64_t load;
#endif

    void SW(const T *A, const T *B, T *C) override;   
};

template<class T>
Winograd<T>::Winograd(i_type n, i_type odd_cap, i_type winograd_cap) : n_(n) {
    bool odd = (n % 2 != 0);
    n += odd;
    bool cap = false;

    while ((n /= 2) != 1) {
        if ((odd && n < odd_cap) || n * 2 < winograd_cap) {
            L_.push_back(std::make_unique<LevelClassic>(n * 2 - odd));
            cap = true;
            break;
        }
        if (odd)
            L_.push_back(std::make_unique<LevelOdd>(n));
        else
            L_.push_back(std::make_unique<LevelEven>(n));

        odd = (n % 2 != 0);
        n += odd;
    }

    if (!cap) {
        L_.push_back(std::make_unique<Level22>());
    }

    for (int i = L_.size() - 1; i > 0; --i) {
        L_[i - 1]->next = &(*L_[i]);
    }

}


template<class T>
void Winograd<T>::Execute(const M &A, const M &B, M &C) {
    if (A.GetCols() != n_ || B.GetCols() != n_ || C.GetCols() != n_ ||
        A.GetRows() != n_ || B.GetRows() != n_ || C.GetRows() != n_) {
        throw std::invalid_argument("Matrix size not match");
    }
    L_[0]->SW(A.Data().data(), B.Data().data(), C.Data().data());
}

template<class T>
void Winograd<T>::Level22::SW(const T *A, const T *B, T *C) {
#ifdef LEVEL_LOAD_TEST__P
    auto time_point = Time::Now();
#endif

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

#ifdef LEVEL_LOAD_TEST__P
    load += Time::Duration<Time::ns>(time_point);
#endif
}

template<class T>
void Winograd<T>::LevelClassic::SW(const T *A, const T *B, T *C) {
#ifdef LEVEL_LOAD_TEST__P
    auto time_point = Time::Now();
#endif
    for (i_type i = 0; i < n; ++i) {
        for (i_type j = 0; j < n; ++j) {
            T sum = 0;
            for (i_type k = 0; k < n; ++k) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
#ifdef LEVEL_LOAD_TEST__P
    load += Time::Duration<Time::ns>(time_point);
#endif
};

template<class T>
void Winograd<T>::LevelEven::SW(const T *A, const T *B, T *C) {

#ifdef LEVEL_LOAD_TEST__P
    auto time_point = Time::Now();
#endif

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

#ifdef LEVEL_LOAD_TEST__P
    level_load[n] += Time::Duration<Time::ns>(time_point);
#endif


    next->SW(A11, B11, R1);
    next->SW(A12, B21, R2);
    next->SW(S4, B22, R3);
    next->SW(A22, T4, R4);
    next->SW(S1, T1, R5);
    next->SW(S2, T2, R6);
    next->SW(S3, T3, R7);

#ifdef LEVEL_LOAD_TEST__P
    time_point = Time::Now();
#endif


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

#ifdef LEVEL_LOAD_TEST__P
    level_load[n] += Time::Duration<Time::ns>(time_point);
#endif

}

template<class T>
void Winograd<T>::LevelOdd::SW(const T *A, const T *B, T *C) {

#ifdef LEVEL_LOAD_TEST__P
    auto time_point = Time::Now();
#endif

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
#ifdef LEVEL_LOAD_TEST__P
    level_load[n] += Time::Duration<Time::ns>(time_point);
#endif

    next->SW(A11, B11, R1);
    next->SW(A12, B21, R2);
    next->SW(S4, B22, R3);
    next->SW(A22, T4, R4);
    next->SW(S1, T1, R5);
    next->SW(S2, T2, R6);
    next->SW(S3, T3, R7);

#ifdef LEVEL_LOAD_TEST__P
    time_point = Time::Now();
#endif


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
#ifdef LEVEL_LOAD_TEST__P
    level_load[n] += Time::Duration<Time::ns>(time_point);
#endif
}

template<class T>
Winograd<T>::PointerBase::PointerBase(i_type n) : n(n) {
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

template<class T>
Winograd<T>::PointerBase::~PointerBase() {
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

} // namespase Basic

} // namespace s21

