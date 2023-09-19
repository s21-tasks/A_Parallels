#include "m_thread.h"

ThreadManager::ThreadManager(int thread_count)
    : ready_(thread_count), thread_count_(thread_count) {
  std::fill(ready_.begin(), ready_.end(), false);
  for (int thc = 0; thc < thread_count; ++thc) {
    threads_.emplace_back([&, thc, thread_count] {
      while (true) {
        {
          std::unique_lock<std::mutex> lock(mtx_);
          cv_.wait(lock, [&] { return (bool)ready_[thc]; });
        }

        if (stopped_) {
          return;
        }

        try {
          func_(thc);
        } catch (std::runtime_error &ex) {
          std::cout << "at " << std::this_thread::get_id() << ' ' << ex.what()
                    << '\n';
        }
        {
          std::unique_lock<std::mutex> lock(mtx_);
          ready_[thc] = false;
          ++finished_;
        }
        cv_.notify_all();
      }
    });
  }
}

void ThreadManager::ShutDown() {
  std::fill(ready_.begin(), ready_.end(), true);
  stopped_ = true;
  cv_.notify_all();
  for (auto &thread : threads_) {
    thread.join();
  }
}

void ThreadManager::Run() {
  std::fill(ready_.begin(), ready_.end(), true);
  cv_.notify_all();
  {
    std::unique_lock<std::mutex> lock(mtx_);
    cv_.wait(lock, [&] { return finished_ == (int)ready_.size(); });
  }
  finished_ = 0;
}

void ThreadManager::Execute(const std::function<void(int)> &func) {
  func_ = func;
  Run();
}

void ThreadManager::LoopExecute(int loop_count,
                                const std::function<void(int)> &func) {
  func_ = LoopThreads(thread_count_, loop_count, func);
  Run();
}

void ThreadManager::DispThreads(int thread_count,
                                const std::function<void(int)> &func) {
  std::vector<std::thread> threads;
  threads.reserve(thread_count);
  for (int thc = 0; thc < thread_count; ++thc) {
    threads.emplace_back(func, thc);
  }
  for (auto &thread : threads) {
    thread.join();
  }
}

std::function<void(int)> ThreadManager::LoopThreads(
    int thread_count, int loop_k, const std::function<void(int)> &func) {
  return [=, &func](int thc) {
    int delta = loop_k / thread_count;
    int end = (thc == thread_count - 1 ? loop_k : (thc + 1) * delta);
    for (int k = thc * delta; k < end; ++k) {
      func(k);
    }
  };
}