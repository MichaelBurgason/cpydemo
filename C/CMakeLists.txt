cmake_minimum_required(VERSION 3.11)

# Include directories
include_directories(include)

# Create static library
add_library(cpydemo_lib SHARED src/CPYDemo.cpp)

# Make include directory visible to other targets
target_include_directories(cpydemo_lib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Set C/C++ standard for the library
set_target_properties(cpydemo_lib PROPERTIES
  C_STANDARD 11
  C_STANDARD_REQUIRED ON
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)

# Add executable
add_executable(CPYDemo src/main.cpp)

# Link the executable with the library
target_link_libraries(CPYDemo PRIVATE cpydemo_lib)

# Set C/C++ standard for the executable
set_target_properties(CPYDemo PROPERTIES
  C_STANDARD 11
  C_STANDARD_REQUIRED ON
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)