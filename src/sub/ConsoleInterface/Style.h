#pragma once

#include <initializer_list>
#include <iostream>
#include <string>

struct Style {
  inline static std::string error_color = "\033[0;31m";
  inline static std::string request_color = "\033[0;35m";
  inline static std::string default_color = "\033[0m";
  inline static std::string exit_word = "exit";
  inline static std::string home_word = "home";

  static void InputRequest(const std::initializer_list<std::string> &requests) {
    std::cout << request_color
              << (requests.size() > 1 ? "Enter separated by a space: "
                                      : "Enter: ");
    auto iter = requests.begin();
    std::cout << *iter;
    while (++iter != requests.end()) {
      std::cout << "; " << *iter;
    }
    std::cout << " or '" << exit_word << " or '" << home_word << "'\n"
              << default_color;
  }

  static void MenuRequest() {
    std::cout << request_color << "Enter number of your choice or '"
              << exit_word << " or '" << home_word << "'\n"
              << default_color;
  }

  static void ErrorPrint(const std::string &error) {
    std::cout << error_color << "Error: \"" << error << "\". Please try again\n"
              << default_color;
  }

  static std::string IncorrectInput(const std::string &input) {
    return "incorrect input " + input;
  }

 private:
};
