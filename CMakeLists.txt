cmake_minimum_required(VERSION 3.8.2)
project(integer_polyomino)
# Tell CMake where to find Python.h
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
# Tell CMake to set the compiler to C++ 17 standards
set (CMAKE_CXX_STANDARD 17)
# Make sure openMP is installed
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set source directory
set (SOURCE_DIR "src/integer_polyomino")
# Set submodule directory
set (GPMAP_DIR "${SOURCE_DIR}/gpmap")
set (ASSEMBLY_DIR "${SOURCE_DIR}/assembly")
# Tell CMake that headers are also in SOURCE_DIR/includes
include_directories(${GPMAP_DIR}/includes)
include_directories(${ASSEMBLY_DIR}/includes)
# Tell CMake where to find other headers library
include_directories(lib/polyomino_core/includes)
include_directories(lib/pybind11/include)

# List the source files
set (SOURCES_GPMAP "${GPMAP_DIR}/mapping.cpp"
            "${GPMAP_DIR}/io.cpp"
            "${GPMAP_DIR}/core_metrics.cpp"
            "${GPMAP_DIR}/metrics.cpp"
            "${GPMAP_DIR}/generate.cpp"
            "${GPMAP_DIR}/duplicate.cpp"
            "${GPMAP_DIR}/gpmap_api.cpp"
            "${ASSEMBLY_DIR}/integer_model.cpp")

set (SOURCES_ASSEMBLY "${ASSEMBLY_DIR}/integer_model.cpp"
            "${ASSEMBLY_DIR}/assembly_api.cpp")


# # Add and compile the test suite for the cpp functions
# SET(TEST_DIR "tests")
# SET(TESTS ${SOURCES}
#     "${TEST_DIR}/test_main.cpp"
#     "${TEST_DIR}/test_assembly.cpp")
#
# # Generate a test executable
# include_directories(lib/Catch2/include)
# add_executable("${PROJECT_NAME}_test" ${TESTS})

# Generate Python module
add_subdirectory(lib/pybind11)
pybind11_add_module(gpmap ${SOURCES_GPMAP} "${GPMAP_DIR}/bindings.cpp")
pybind11_add_module(assembly ${SOURCES_ASSEMBLY} "${ASSEMBLY_DIR}/bindings.cpp")
