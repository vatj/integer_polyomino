cmake_minimum_required(VERSION 3.8.2)
project(integer_polyomino)
# Tell CMake where to find Python.h
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
# Tell CMake to set the compiler to C++ 17 standards
SET(CMAKE_CXX_STANDARD 17)
# Make sure openMP is installed
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set source directory
SET(SOURCE_DIR "src/integer_polyomino")
# Tell CMake that headers are also in SOURCE_DIR/includes
include_directories(${SOURCE_DIR}/includes)
# Tell CMake where to find other headers library
include_directories(lib/polyomino_core/includes)
include_directories(lib/pybind11/include)
# List the source files
SET(SOURCES "${SOURCE_DIR}/integer_model.cpp"
            "${SOURCE_DIR}/genotype_phenotype.cpp"
            "${SOURCE_DIR}/genotype_iofunc.cpp"
            "${SOURCE_DIR}/genotype_core_metrics.cpp"
            "${SOURCE_DIR}/genotype_metrics.cpp"
            "${SOURCE_DIR}/genotype_generate.cpp"
            "${SOURCE_DIR}/genotype_duplicate.cpp"
            "${SOURCE_DIR}/genotype_api.cpp")


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
pybind11_add_module(integer_polyomino ${SOURCES} "${SOURCE_DIR}/bindings.cpp")
