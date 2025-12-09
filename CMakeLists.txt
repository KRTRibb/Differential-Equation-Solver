cmake_minimum_required(VERSION 3.20)
project(diffeq_demo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

set(PYTHON_VENV_PATH "${CMAKE_SOURCE_DIR}/.venv")
set(Python3_ROOT_DIR "${PYTHON_VENV_PATH}")
set(Python3_FIND_VIRTUALENV ONLY)
set(Python3_FIND_STRATEGY LOCATION)

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

message(STATUS "🔹 Using Python: ${Python3_EXECUTABLE}")
message(STATUS "🔹 Python include dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "🔹 NumPy include dir: ${Python3_NumPy_INCLUDE_DIRS}")

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED
    HINTS /opt/homebrew/opt/pybind11/share/cmake/pybind11
          /usr/local/opt/pybind11/share/cmake/pybind11
)

include_directories(${pybind11_INCLUDE_DIRS})

add_executable(diffeq_demo
    src/main.cpp
    src/solver.cpp
    src/Steppers/euler_stepper.cpp
)

target_include_directories(diffeq_demo PRIVATE include)

pybind11_add_module(diffeqpy
    src/bindings.cpp
    src/solver.cpp
    src/Steppers/euler_stepper.cpp
    src/Steppers/midpoint_stepper.cpp
    src/Steppers/rk3_stepper.cpp
    src/Steppers/rk4_stepper.cpp
    src/Steppers/rk5_stepper.cpp
    src/pde.cpp
)

target_include_directories(diffeqpy PRIVATE
    include
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
)

target_link_libraries(diffeqpy PRIVATE Python3::Python)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
        "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE VENV_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "🔹 Installing diffeqpy to: ${VENV_SITE_PACKAGES}")

if(EXISTS "${VENV_SITE_PACKAGES}")
    add_custom_command(TARGET diffeqpy POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
            "$<TARGET_FILE:diffeqpy>"
            "${VENV_SITE_PACKAGES}/"
        COMMENT "📦 Installing diffeqpy into active Python venv site-packages"
    )
else()
    message(WARNING "⚠️ Could not detect site-packages directory at ${VENV_SITE_PACKAGES}. Skipping install step.")
endif()

add_custom_target(count
    COMMAND find ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include -name "*.cpp" -o -name "*.hpp" -o -name "*.h" |
            xargs wc -l
    COMMENT "Counting total lines of C++ code in src/ and include/"
)
