cmake_minimum_required(VERSION 3.16)
project(SupCompMod24)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMP package
find_package(OpenMP REQUIRED)


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
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-Werror")
enable_cxx_compiler_flag_if_supported("-Wpedantic")
enable_cxx_compiler_flag_if_supported("-Wnon-virtual-dtor")

add_executable(linear src/main.cpp)
add_executable(openmp src/mainOpenMP.cpp)

# Link OpenMP to the executable
if(OpenMP_CXX_FOUND)
    target_compile_options(openmp PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(openmp PUBLIC OpenMP::OpenMP_CXX)
endif()