cmake_minimum_required(VERSION 3.12)
project(Tutorial LANGUAGES CXX)

if(UNIX AND NOT APPLE)
    set(CMAKE_C_COMPILER "/usr/bin/gcc")
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
elseif(APPLE)
    set(CMAKE_C_COMPILER "clang")
    set(CMAKE_CXX_COMPILER "clang++")
else()
    message(WARNING "Only support LINUX and APPLE for now.")
endif()

# for color output only
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

SET(CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,${CMAKE_SOURCE_DIR}/build")

# find the following packages
find_package(Eigen3 REQUIRED HINTS ${Eigen3_HINT})
# for using casadi
find_package(casadi REQUIRED)
# For multi threading
find_package(Threads REQUIRED)

# find_package(Ipopt REQUIRED)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# set the g++ compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")


if(CMAKE_CXX_COMPILER_LOADED)
    message(STATUS "Compiler Path: ${CMAKE_CXX_COMPILER}")
    message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")
endif()

# set the Ipopt library and include directory
# set(LIBRARY_DIRS /usr/local/lib)
# if(UNIX AND NOT APPLE)
#     set(IPOPT_INCLUDE_DIRS /usr/include/coin-or)
# elseif(APPLE)
#     set(IPOPT_INCLUDE_DIRS /usr/local/include/coin-or)
# endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(IPOPT REQUIRED ipopt)
# pkg_check_modules(SNOPT REQUIRED libsnopt7_cpp)
# include_directories(${IPOPT_INCLUDE_DIRS})
# link_directories(${IPOPT_LIBRARY_DIRS})

# set casadi include directory
set(CASADI_INCLUDE_DIR /usr/local/include/casadi)

# find casadi library
find_library(CASADI_LIBRARY
    NAMES casadi
    HINTS ${CASADI_INCLUDE_DIR}/../lib $ENV{CASADI_PREFIX}/lib)
if(CASADI_LIBRARY)
    set(CASADI_LIBRARIES ${CASADI_LIBRARIES} ${CASADI_LIBRARY})
endif()



# include header directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    SYSTEM ${IPOPT_INCLUDE_DIRS}
    SYSTEM ${CASADI_INCLUDE_DIR}
    SYSTEM ${SNOPT_INCLUDE_DIRS}
    ${CODE_GEN_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIRS})
include_directories(/home/basit/workspace/libsnopt7_cpp)


# library directories
link_directories(${LIBRARY_DIRS}
    ${IPOPT_LIBRARY_DIRS})
link_directories(/home/basit/workspace/libsnopt7_cpp)


# find all the header files
file(GLOB HEADER_FILES_HPP ${CMAKE_SOURCE_DIR}/include/*.hpp)
file(GLOB HEADER_FILES_H ${CMAKE_SOURCE_DIR}/include/*.h)


# add executable files for test

# casadi + ipopt test
# add_executable(example_racecar example/example_racecar.cpp)
# target_link_libraries(example_racecar ${CASADI_LIBRARIES} ipopt)

# # casadi basic operations
# add_executable(example_basic_operations example/example_basic_operations.cpp)
# target_link_libraries(example_basic_operations ${CASADI_LIBRARIES})

# # casadi function C code auto generation
# add_executable(example_code_gen example/example_code_gen.cpp)
# target_link_libraries(example_code_gen ${CASADI_LIBRARIES})

# a test file to use auto generated code library
# compile the test file


# # solve a NLP given a compiled shared library
# add_executable(example_use_nlp_external example/example_use_nlp_external.cpp)
# target_link_libraries(example_use_nlp_external ${CASADI_LIBRARIES})

# Define the library target
# add_library(casadiprint SHARED src/casadiprint.cpp)

# Link the CasADi library to your target
# target_link_libraries(casadiprint PRIVATE casadi)

# Pass the compiler flags necessary for CasADi
# set_target_properties(casadiprint PROPERTIES
#                       COMPILE_FLAGS "${CASADI_CFLAGS_OTHER}")

add_executable(usecodegen src/usecodegen.cpp)
set_target_properties(usecodegen PROPERTIES INSTALL_RPATH "/home/basit/workspace/libsnopt7_cpp")
target_link_libraries(usecodegen ${CASADI_LIBRARIES} ipopt -Llibsnopt7_cpp)

# compile a nonlinear programming to a shared library
add_executable(nmpccodegen src/nmpccodegen.cpp) 
set_target_properties(nmpccodegen PROPERTIES INSTALL_RPATH "/home/basit/workspace/libsnopt7_cpp")
target_link_libraries(nmpccodegen ${CASADI_LIBRARIES} ipopt -Llibsnopt7_cpp)