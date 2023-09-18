#pragma once

#include "winograd.h"
#include "../sub/thread_pipeline/thread_pipeline.h"

namespace s21 {


namespace PipeLine {

template<class T>
class Winograd {
    using M = Matrix<T>;
    using i_type = typename M::i_type;
    public:
        Winograd(const std::vector<M> &workers, i_type N);
        void Execute(M &m);


        static void Mul(const M &A, const M &B, M &C);
    private:
        Pipeline<M> pipeline_{};
        std::vector<M> workers_;
};

template<class T>
Winograd<T>::Winograd(const std::vector<M> &workers, i_type N) :
        workers_(workers) {
    
    i_type s = workers.size();

    for (i_type k = 0; k < s; ++k) {
        std::shared_ptr<Basic::Winograd<T>> W_ptr = std::make_shared<Basic::Winograd<T>>(N);
        pipeline_.AddStage([W_ptr, this, k](M &m) {
            W_ptr->Execute(m, workers_[k], m);
        });
        std::cout << "Stage added\n";
    }
}

template<class T>
void Winograd<T>::Execute(M &m) {
    pipeline_.Process(m);
}

template<class T>
void Winograd<T>::Mul(const M &A, const M &B, M &C) {
    Basic::Winograd<T>::Mul(A, B, C);
}

} // namespace PipeLine

} // namespace s21
