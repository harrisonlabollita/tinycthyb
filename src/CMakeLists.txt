# Start configuration
cmake_minimum_required(VERSION 2.8.12)
project(tinycthyb CXX)
set(CMAKE_BUILD_TYPE Release)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# Create executable
add_executable(main main.cpp)

# Linking and include info
target_link_libraries(main triqs)
triqs_set_rpath_for_target(main)
