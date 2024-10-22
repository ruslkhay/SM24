cmake_minimum_required(VERSION 3.17)
project(SupCompMod24)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMP package
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Specify the directories for include files
include_directories(src)
# Collect all source files from the src directory
file(GLOB_RECURSE SRC_FILES src/*.cpp src/*.hpp)

include(CheckCXXCompilerFlag)
function(enable_cxx_compiler_flag_if_supported flag)
  string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
  if(flag_already_set EQUAL -1)
    check_cxx_compiler_flag("${flag}" flag_supported)
    if(flag_supported)
      set(CMAKE_CXX_FLAGS
          "${CMAKE_CXX_FLAGS} ${flag}"
          PARENT_SCOPE)
    endif()
    unset(flag_supported CACHE)
  endif()
endfunction()
enable_cxx_compiler_flag_if_supported("-O2")
enable_cxx_compiler_flag_if_supported("-g")
enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Werror")

# Google benchmark
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.0)
FetchContent_MakeAvailable(googletest)

include(FetchContent)
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.8.4)
FetchContent_MakeAvailable(googlebenchmark)

add_executable(bench benchmarks.cpp)
target_sources(bench PRIVATE ${SRC_FILES} )
target_link_libraries(bench PRIVATE benchmark::benchmark OpenMP::OpenMP_CXX)
