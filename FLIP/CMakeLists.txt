cmake_minimum_required(VERSION 3.10)
project(FLIP)

if(FALSE)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
endif()

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

target_link_libraries(run ${OpenCV_LIBS} fluid)