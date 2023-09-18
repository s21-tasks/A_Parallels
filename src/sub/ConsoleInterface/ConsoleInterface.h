#pragma once

#include <functional>
#include <initializer_list>
#include <limits>
#include <list>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Style.h"

class Menu;
class Option;
class ConsoleInterface;

class AbsMenu {
 protected:
  std::string name_;
  ConsoleInterface *ci_;

 public:
  virtual void Action() = 0;
  AbsMenu(const std::string &name, ConsoleInterface *ci);
  std::string &GetName();
  virtual ~AbsMenu() = default;
  void GoHome();
};

class Input : public AbsMenu {
 public:
  Input(ConsoleInterface *ci, AbsMenu *next_menu, const std::string &name = "");
  virtual ~Input() = default;

 protected:
  void Action() override;
  virtual void Func() = 0;

  template <class T>
  T Read();
  bool Allowed();

 private:
  bool exit_ = false, home_ = false;
  AbsMenu *next_menu_;
};

class Menu : public AbsMenu {
 public:
  Menu(const std::initializer_list<std::string> &options,
       const std::string &name, ConsoleInterface *ci);

  void Action() override;
  void Connect(
      int k, AbsMenu *next_menu,
      const std::function<void(void)> &action = [] {});

 private:
  class Option {
    friend class Menu;

   public:
    Option() {}
    Option(const std::string &name);

    void operator()();

   private:
    std::string name_;
    std::function<void(void)> action_;
    AbsMenu *next_menu_;
  };

  std::vector<Option> options_;
  const std::string name_;
};

class ConsoleInterface {
 public:
  ConsoleInterface() {}

  Menu *AddMenu(const std::initializer_list<std::string> &options,
                const std::string &name = "");
  AbsMenu *AddInput(Input *input);
  void Start();
  virtual ~ConsoleInterface() = default;

 private:
  friend AbsMenu;
  std::list<std::unique_ptr<AbsMenu>> menus_;
};

template <class T>
T Input::Read() {
  if (Allowed()) {
    std::string word;
    std::cin >> word;
    T value;
    if (word == Style::exit_word) {
      exit_ = true;
    } else if (word == Style::home_word) {
      home_ = true;
    } else {
      std::istringstream iss(word);
      if (!(iss >> value)) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        throw std::runtime_error(Style::IncorrectInput(word));
      }
      return value;
    }
  }
  return T();
}

template <class T>
class OneValueInput : public Input {
 public:
  OneValueInput(ConsoleInterface *ci, AbsMenu *next_menu,
                std::function<void(T &)> func, const std::string request,
                const std::string &name = "")
      : Input(ci, next_menu, name), request_(request), func_(func) {}

  void Func() override {
    Style::InputRequest({request_});
    T value = Read<T>();
    if (Allowed()) {
      func_(value);
    }
  }

 private:
  std::string request_;
  std::function<void(T &)> func_;
};
