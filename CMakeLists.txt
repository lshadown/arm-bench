cmake_minimum_required(VERSION 3.28)
project(arm-bench)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -funroll-loops")
find_package(OpenMP REQUIRED)

add_executable(test
        test.cpp
)

if (OpenMP_FOUND)
    target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)
    if (APPLE)
        target_compile_options(test PUBLIC -Xpreprocessor -fopenmp)
        target_link_libraries(test PUBLIC /opt/homebrew/opt/libomp/lib/libomp.dylib)
        target_include_directories(test PUBLIC /opt/homebrew/opt/libomp/include)
    endif()
else()
    message(FATAL_ERROR "OpenMP not found")
endif()

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")