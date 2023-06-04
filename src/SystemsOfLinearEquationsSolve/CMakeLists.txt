cmake_minimum_required(VERSION 3.22)
project(gauss VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#add_subdirectory(GraphAlgorithms ${CMAKE_BINARY_DIR}/GraphAlgorithms)
#add_subdirectory(Submodules/ConsoleInterface)


add_executable( ${PROJECT_NAME}
        SLAEGauss.cc
        SLAEGauss.h
        )


#target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/TslCompare)
#target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/Interface)

#target_link_libraries(${PROJECT_NAME} PRIVATE s21_graph_algorithms.a)
#target_link_libraries(${PROJECT_NAME} PRIVATE console_interface)
