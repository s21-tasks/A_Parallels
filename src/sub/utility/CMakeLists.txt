cmake_minimum_required(VERSION 3.0)

project(utility VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_library(${PROJECT_NAME} STATIC
        utility.h
        m_sstr.h
        m_time.h
        dummy.cc
        m_random.h)

set_target_properties(${PROJECT_NAME} PROPERTIES
OUTPUT_NAME cpp_utilites
PREFIX ""
)


target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
