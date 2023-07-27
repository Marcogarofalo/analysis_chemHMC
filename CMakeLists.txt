cmake_minimum_required(VERSION 3.13)

project(analysis_chemHMC CXX C)


IF(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -static-libasan -g   ")#-fsanitize=float-divide-by-zero
  set(LDFLAGS "${LDFLAGS}  -fsanitize=address")
  set(CMAKE_LDFLAGS_FLAGS  " ${CMAKE_LDFLAGS_FLAGS} ${LDFLAGS} ")
ENDIF()

#############################################
# load analysis_program
#############################################
find_package(analysis_program)
include_directories(${ANALYSIS_PROGRAM_INCLUDE_DIR})
MESSAGE("ANALYSIS_PROGRAM_INCLUDE_DIR   ${ANALYSIS_PROGRAM_INCLUDE_DIR}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARIES: ${ANALYSIS_PROGRAM_LIBRARIES}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARY_DIR: ${ANALYSIS_PROGRAM_LIBRARY_DIR}")
IF(NOT analysis_program_FOUND)
    MESSAGE("find_package faild to find yaml-cpp: you need to do it manually")
    MESSAGE("hint:  analysis_program_FOUND_DIR or CMAKE_PREFIX_PATH need to be global path or relative to the source dir")
    MESSAGE(FATAL_ERROR)
ENDIF()
find_package(analysis_program)

find_package(Eigen3 REQUIRED) # we need eigen

#############################################
# add a local library for this project
#############################################
add_library(
    functions_analysis_chemHMC STATIC
    functions_analysis_chemHMC.hpp
    functions_analysis_chemHMC.cpp
)
target_link_libraries(functions_analysis_chemHMC PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
target_include_directories(functions_analysis_chemHMC PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(analysis_chemHMC analysis_chemHMC.cpp)
target_link_libraries(analysis_chemHMC PUBLIC functions_analysis_chemHMC)
target_link_libraries(analysis_chemHMC PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
