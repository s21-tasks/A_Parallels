cmake_minimum_required(VERSION 3.0)

project(console_interface VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_library(${PROJECT_NAME} STATIC
        ConsoleInterface.h
        Style.h
        ConsoleInterface.cc)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})