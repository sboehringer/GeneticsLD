cmake_minimum_required(VERSION 2.6)
project(haplotype_cpp)

set(CMAKE_CXX_FLAGS "-std=c++11")
include_directories(/usr/local/lib64/R-3.0.1/library/Rcpp/include /usr/local/lib64/R-3.0.1/include/)
link_directories(/usr/local/lib64/R-3.0.1/lib /usr/local/lib64/R-3.0.1/library/Rcpp/lib)

add_executable(haplotypes_cpp main.cpp parameter.cpp parameterhelper.cpp parametrizer.cpp)
add_library(haplotypes_rcpp SHARED parameter.cpp parameterhelper.cpp parametrizer.cpp parametrizerrcpp.cpp)
target_link_libraries(haplotypes_rcpp R Rcpp)

install(TARGETS haplotypes_cpp RUNTIME DESTINATION bin)
