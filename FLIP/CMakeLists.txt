cmake_minimum_required(VERSION 3.10)
project(FLIP)

if(FALSE)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
endif()

find_package(OpenMP)

add_executable(run run.cpp)

add_library(fluid)

target_compile_features(fluid
	PUBLIC cxx_std_17)

target_include_directories(fluid
	PUBLIC include/)

target_sources(fluid
	PRIVATE
		"src/fluid.cpp"
		"src/grid.cpp"
		"src/pressure_solver.cpp")

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(run PUBLIC OpenMP::OpenMP_CXX fluid)
  target_compile_options(run PRIVATE "${OpenMP_CXX_FLAGS}")
else()
  message(WARNING "OpenMP not found. Building without OpenMP support.")
  target_link_libraries(run ${OpenCV_LIBS} fluid)
endif()