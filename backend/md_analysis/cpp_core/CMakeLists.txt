cmake_minimum_required(VERSION 3.12)
project(molecular_interactions)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find pybind11
find_package(pybind11 CONFIG REQUIRED)

# Create the pybind11 module
pybind11_add_module(molecular_interactions 
    pybind11_module.cpp
    radii_lists.cpp
    nonbond_interactions.cpp
    bonded_interactions.cpp
)

# Set optimization flags
if(MSVC)
    target_compile_options(molecular_interactions PRIVATE /O2 /W4)
else()
    target_compile_options(molecular_interactions PRIVATE -O3 -Wall -Wextra)
endif()

# Installation
install(TARGETS molecular_interactions LIBRARY DESTINATION .)
