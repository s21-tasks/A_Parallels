#pragma once

#include <filesystem>
#include <iostream>
#include <iterator>
#include <sstream>
#include <type_traits>

namespace s21 {

struct SStr {
  template <const char Delimiter, class... Args>
  static std::string Fill(const Args &...args) {
    std::ostringstream sstr;
    int i = 0;
    ((sstr << std::dec << args << (++i == sizeof...(args) ? '\0' : Delimiter)),
     ...);
    return sstr.str();
  }

  template <class... Args>
  static std::string Fill(const Args &...args) {
    std::ostringstream sstr;
    ((sstr << std::dec << args), ...);
    return sstr.str();
  }

  template <const char Delimiter = '\n', class... Args>
  static void Print(const Args &...args) {
    int i = 0;
    (
        [&] {
          PrintP(args);
          PrintStyle::stream << (++i == sizeof...(args) ? '\0' : Delimiter);
        }(),
        ...);
    PrintStyle::stream << '\n';
  }

  template <class T>
  static void Print(const T &arg) {
    PrintP(arg);
    PrintStyle::stream << '\n';
  }

  struct PrintStyle {
    inline static std::string container_delimiter = ", ";
    inline static std::string pair_delimiter = ", ";
    inline static std::ostream &stream = std::cout;

    template <class T>
    static typename std::enable_if<std::is_same<T, std::string>::value,
                                   std::string>::type
    Wrapper(const T &l) {
      return "\"" + l + "\"";
    }

    template <class T>
    static typename std::enable_if<!std::is_same<T, std::string>::value,
                                   const T &>::type
    Wrapper(const T &l) {
      return l;
    }
  };

 private:
  template <typename T>
  struct is_container {
    template <typename U>
    static constexpr bool is_container_helper(
        U *, decltype(std::declval<U>().begin()) * = nullptr,
        decltype(std::declval<U>().end()) * = nullptr) {
      return true;
    }

    template <typename>
    static constexpr bool is_container_helper(...) {
      return false;
    }

    static constexpr bool value = is_container_helper<T>(nullptr);
  };

  template <class Container>
  static typename std::enable_if<is_container<Container>::value, void>::type
  PrintP(const Container &c) {
    auto i = c.begin();
    if (i == c.end()) {
      return;
    }
    PrintStyle::stream << '{';
    while (i != --c.end()) {
      PrintP(*(i++));
      PrintStyle::stream << PrintStyle::container_delimiter;
    }
    PrintP(*(i++));
    PrintStyle::stream << "}";
  }

  template <class T>
  static typename std::enable_if<!is_container<T>::value, void>::type PrintP(
      const T &data) {
    PrintStyle::stream << data;
  }

  template <class A, class B>
  static void PrintP(const std::pair<A, B> &p) {
    // PrintStyle::stream << '[' << p.first << PrintStyle::pair_delimiter <<
    // p.second << "]";
    PrintStyle::stream << '[';
    PrintP(p.first);
    PrintStyle::stream << PrintStyle::pair_delimiter;
    PrintP(p.second);
    PrintStyle::stream << "]";
  }

  template <class A>
  static void PrintP(const std::pair<A, A> &p) {
    // PrintStyle::stream << '[' << p.first << PrintStyle::pair_delimiter <<
    // p.second << "]";
    PrintStyle::stream << '[';
    PrintP(p.first);
    PrintStyle::stream << PrintStyle::pair_delimiter;
    PrintP(p.second);
    PrintStyle::stream << "]";
  }

 public:
  static std::string RelativePath(const char *FILE,
                                  const std::string &path_from_source) {
    using namespace std::filesystem;
    return relative(path(FILE).parent_path(), current_path()).string() +
           path_from_source;
  }
};

}  // namespace s21
