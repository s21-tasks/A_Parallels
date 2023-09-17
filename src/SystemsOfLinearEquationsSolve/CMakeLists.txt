cmake_minimum_required(VERSION 3.22)
project(gauss VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../sub/matrix ${CMAKE_BINARY_DIR}/matrix)
add_subdirectory(../sub/utility ${CMAKE_BINARY_DIR}/utility)
add_subdirectory(../sub/ConsoleInterface ${CMAKE_BINARY_DIR}/interface)


set(SOURCES
        SLAEGauss.cc
        SLAEGauss.h
        tests/GaussSpeedTest.cc
        tests/GaussSpeedTest.h
        tests/main.cc
    )

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC utility)
target_link_libraries(${PROJECT_NAME} PUBLIC matrix)
target_link_libraries(${PROJECT_NAME} PRIVATE console_interface)
