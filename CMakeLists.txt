cmake_minimum_required( VERSION 2.8 )

project( calcPi )

set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules" )
set( MAIN calcPi.cpp)

find_package(GMP REQUIRED)

## Set Compiler Flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -x c++ -pie")
add_definitions("-Wall -x c++")

add_executable(calcPi.exec ${MAIN})
target_link_libraries(calcPi.exec gmp gmpxx)