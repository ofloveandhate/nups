cmake_minimum_required (VERSION 2.6)
project (nups)

# The version number.
set (nups_VERSION_MAJOR 0)
set (nups_VERSION_MINOR 1)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/nups_config.h.in"
#  "${PROJECT_BINARY_DIR}/nups_config.h"
#  )

include_directories (include)

set(NUPS_HEADERS
	include/nups/nups.hpp
	include/nups/factor.hpp
	include/nups/linear_solve.hpp
	include/nups/polynomial_solve.hpp
	include/nups/predict.hpp
	include/nups/type_traits.hpp
	)


set(NUPS_SOURCES
	src/nups.cpp
	)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
	
include_directories("${PROJECT_BINARY_DIR}")

add_library(nups ${NUPS_HEADERS} ${NUPS_SOURCES})



IF (BUILD_TESTING)

	#Setup CMake to run tests
	enable_testing()

	#Prep for compiling against boost
	find_package(Boost COMPONENTS unit_test_framework REQUIRED)
	include_directories (${Boost_INCLUDE_DIRS})

	add_executable(nups_test test/nups_test.cpp) 
	target_link_libraries(nups_test ${Boost_LIBRARIES} nups)
	add_test(NUPS nups_test)

ENDIF(BUILD_TESTING)







