#pragma once

#include <iostream>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>
#include <vector>

class ThreadManager {
    public:
        ThreadManager(int thread_count);
        void ShutDown();
        void Execute(const std::function<void(int)> &func);
        void LoopExecute(int loop_count, const std::function<void(int)> &func);
        ~ThreadManager() { ShutDown(); }

        void SetFunction(const std::function<void(int)> &func);
        void SetLoopFunction(int loop_count, const std::function<void(int)> &func);

        static std::function<void(int)> LoopThreads(int thread_count, int loop_k, const std::function<void(int)> &func);
        static void DispThreads(int thread_count, const std::function<void(int)> &func);

        void Run();
    
    private:
        std::vector<std::thread> threads_;
        std::mutex mtx_;
        std::condition_variable cv_;
        std::atomic_int finished_ = 0;
        std::vector<std::atomic_bool> ready_;
        std::atomic_bool stopped_ = false;
        std::function<void(int)> func_;
        int thread_count_;

        // void Run();
};
