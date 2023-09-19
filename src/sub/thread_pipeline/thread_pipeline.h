#pragma once

#include <any>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

template <typename T>
class Pipeline {
  using func_type = std::function<void(T&)>;
  using check_type = std::function<bool(T&)>;

  class PipelineElement {
   public:
    // Удаляем конструктор по умолчанию
    PipelineElement() = delete;

    // Конструктор копирования
    PipelineElement(PipelineElement const& other) {
      func = other.func;
      check = other.check;
      t = new std::thread(&PipelineElement::Worker, this);
    }

    // Оператор присваивания копированием
    PipelineElement& operator=(PipelineElement const& other) {
      if (this == &other) return *this;
      this = new PipelineElement(other);
      return *this;
    }

    // Запрещаем перемещение объектов PipelineElement
    PipelineElement& operator=(PipelineElement&& other) = delete;
    PipelineElement(PipelineElement&& other) = delete;

    // Конструктор класса PipelineElement
    PipelineElement(const func_type f, const check_type ch)
        : func(f), check(ch) {
      t = new std::thread(&PipelineElement::Worker, this);
    }

    // Деструктор класса PipelineElement
    ~PipelineElement() {
      Finish();
      delete next;
    }

    // Обработка данных в элементе конвейера
    void Process(T& data) {
      static std::atomic_int i = 0;
      if (!check || check(data)) {
        {
          std::unique_lock<std::mutex> lock(mtx);
          task_queue.push(data);
        }
        cv.notify_one();
      } else if (next) {
        next->Process(data);
      }
    }

    // Установка следующего элемента конвейера
    void SetNext(PipelineElement* elem) { next = elem; }

    // Получение указателя на следующий элемент конвейера
    inline PipelineElement* GetNext() { return next; }

   private:
    // Рабочая функция элемента конвейера
    void Worker() {
      while (true) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return is_finished || !task_queue.empty(); });
        if (is_finished && task_queue.empty()) {
          return;
        }
        T& data = task_queue.front();
        task_queue.pop();
        lock.unlock();
        try {
          if (!func) {
            throw std::runtime_error("func is nullptr");
          }
          func(data);
        } catch (std::exception& e) {
          std::cerr << e.what() << std::endl;
          continue;
        }

        if (next != nullptr) {
          next->Process(data);
        }
      }
    }

    // Завершение работы элемента конвейера
    void Finish() {
      {
        std::unique_lock<std::mutex> lock(mtx);
        is_finished = true;
      }
      cv.notify_all();
      t->join();
      // next->Finish();
      delete t;
    }

    func_type func;
    PipelineElement* next;
    check_type check;
    std::queue<std::reference_wrapper<T>> task_queue;
    std::condition_variable cv;
    std::mutex mtx;
    bool is_finished = false;
    std::thread* t;
  };

 public:
  // Конструктор класса Pipeline
  Pipeline() = default;

  // Запрещаем копирование и перемещение объектов Pipeline
  Pipeline(Pipeline& other) = delete;
  Pipeline(Pipeline&& other) = delete;
  Pipeline& operator=(Pipeline& other) = delete;
  Pipeline& operator=(Pipeline&& other) = delete;

  // Деструктор класса Pipeline
  ~Pipeline() {
    if (first) {
      delete first;
    }
  }

  // Добавление нового этапа в конвейер
  void AddStage(const func_type func_, const check_type check = nullptr) {
    auto* new_element = new PipelineElement(func_, check);
    if (!first) {
      first = new_element;
      last = new_element;
    } else {
      last->SetNext(new_element);
      last = new_element;
    }
  }

  // Обработка данных в конвейере
  void Process(T& data) {
    if (first) first->Process(data);
  }

 private:
  PipelineElement* first;
  PipelineElement* last;
};