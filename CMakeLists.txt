cmake_minimum_required(VERSION 3.11)
# Force C++ compiler to be set to a specific path
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
project(CPYDemo)

# Set path to external dependencies
set(EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)

# Find Python development libraries
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Add pybind11 subdirectory
set(PYBIND11_PYTHON_VERSION 3)
add_subdirectory(${EXTERNAL_DIR}/pybind11)

# Add the C++ component subdirectory
add_subdirectory(C)

# Create Python module with new name
pybind11_add_module(cpyxdemo_python python/cpyxdemo_module.cpp)

# Link against the CPYDemo library
target_link_libraries(cpyxdemo_python PRIVATE cpydemo_lib)

# Set C++ standard for the Python module
set_target_properties(cpyxdemo_python PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)

# Rename the output Python module to cpyxdemo
set_target_properties(cpyxdemo_python PROPERTIES OUTPUT_NAME "cpyxdemo")