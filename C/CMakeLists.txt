cmake_minimum_required(VERSION 3.11)

# Include directories
include_directories(include)

# Find Python3
find_package(Python3 COMPONENTS Development REQUIRED)

# Find OpenSimplexNoise source files
file(GLOB OPENSIMPLEX_SOURCES 
  "${CMAKE_SOURCE_DIR}/external/openSimplexNoise/OpenSimplexNoise/*.cpp"
)

# Create shared library
add_library(cpydemo_lib SHARED 
            src/CPYDemo.cpp
            #src/TerrainGen.cpp
            src/TerrainGen2DAnimation.cpp
            ${OPENSIMPLEX_SOURCES}
            )

# Make include directory visible to other targets
target_include_directories(cpydemo_lib PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/external
  ${CMAKE_SOURCE_DIR}/external/pybind11/include
  ${Python3_INCLUDE_DIRS}
)

# Link against Python libraries
target_link_libraries(cpydemo_lib PRIVATE
  ${Python3_LIBRARIES}
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
target_link_libraries(CPYDemo PRIVATE 
  cpydemo_lib
  ${Python3_LIBRARIES}
)

# Set C/C++ standard for the executable
set_target_properties(CPYDemo PROPERTIES
  C_STANDARD 11
  C_STANDARD_REQUIRED ON
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)