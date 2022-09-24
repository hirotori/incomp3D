cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

project(Incomp3D Fortran)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

add_subdirectory(lib)

# システムへのインストールは行わない.
# src/メイン関数はプロジェクトに含めない.