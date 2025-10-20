cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0148 NEW)
project(diffeq_demo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---- main demo executable ----
add_executable(diffeq_demo
    src/main.cpp
    src/Steppers/euler_stepper.cpp
    src/solver.cpp
)

target_include_directories(diffeq_demo PRIVATE include)

add_custom_target(count
    COMMAND find ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include -name "*.cpp" -o -name "*.hpp" -o -name "*.h" |
            xargs wc -l
    COMMENT "Counting total lines of C++ code in src/ and include/"
)

# ---- pybind11 bindings ----
# Try to find pybind11 (Homebrew installs it in /opt/homebrew/opt/pybind11)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED
    HINTS /opt/homebrew/opt/pybind11/share/cmake/pybind11
          /usr/local/opt/pybind11/share/cmake/pybind11
)

include_directories(${pybind11_INCLUDE_DIRS})

# ---- Find Python and NumPy ----
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

message(STATUS "Using Python3: ${Python3_EXECUTABLE}")
message(STATUS "Python include dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "NumPy include dir: ${Python3_NumPy_INCLUDE_DIRS}")

# ---- Python module target ----
pybind11_add_module(diffeqpy
    src/bindings.cpp
    src/solver.cpp
    src/Steppers/euler_stepper.cpp
    src/Steppers/midpoint_stepper.cpp
    src/Steppers/rk3_stepper.cpp
)

target_include_directories(diffeqpy PRIVATE
    include
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
)

# ---- Output to 'python/' folder ----
set_target_properties(diffeqpy PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python
)

# Optional: if you build multiple configs (Debug/Release)
foreach(CONFIG_TYPE IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
    set_target_properties(diffeqpy PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY_${CONFIG_TYPE} ${CMAKE_SOURCE_DIR}/python
    )
endforeach()

# ---- Link Python ----
target_link_libraries(diffeqpy PRIVATE Python3::Python)
