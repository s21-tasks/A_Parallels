cmake_minimum_required(VERSION 3.0)

project(matrix VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(${PROJECT_NAME} STATIC
        matrix.h
        dummy.cc)

set_target_properties(${PROJECT_NAME} PROPERTIES
OUTPUT_NAME matrix
PREFIX ""
)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})