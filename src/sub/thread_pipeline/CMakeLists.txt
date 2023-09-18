cmake_minimum_required(VERSION 3.0)

project(thread_pipeline VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_compile_options(-Wall -Wextra -Werror)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
        thread_pipeline.h
        dummy.cc
        )

add_library(${PROJECT_NAME} STATIC ${SOURCES})

add_executable(thread_pipeline_test
        ${SOURCES}
        test.cc)


set_target_properties(${PROJECT_NAME} PROPERTIES
        OUTPUT_NAME thread_pipeline
        PREFIX ""
        )

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
