cmake_minimum_required(VERSION 3.17)
project(SupCompMod24)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMP package
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(MPI REQUIRED) #make it REQUIRED, if you want
# include_directories(SYSTEM ${MPI_INCLUDE_PATH})

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

add_executable(main main.cpp)
target_sources(main PRIVATE ${SRC_FILES} )
target_link_libraries(main PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(main PRIVATE "${OpenMP_CXX_FLAGS}")

add_executable(mpi mainMPI.cpp)
target_sources(mpi PRIVATE ${SRC_FILES} )
target_link_libraries(mpi PUBLIC MPI::MPI_CXX)
# target_link_libraries(mpi PRIVATE ${MPI_CXX_LIBRARIES})