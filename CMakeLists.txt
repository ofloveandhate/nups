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

IF (BUILD_SHARED_LIBS)
	add_library(nups SHARED ${NUPS_HEADERS} ${NUPS_SOURCES})
ELSE ()
	add_library(nups STATIC ${NUPS_HEADERS} ${NUPS_SOURCES})
ENDIF(BUILD_SHARED_LIBS)


install(
		TARGETS nups 
		ARCHIVE DESTINATION "lib"
		LIBRARY DESTINATION "lib"
		COMPONENT library)

install(
		FILES ${NUPS_HEADERS} 
		DESTINATION "include/nups")


enable_testing()
IF (BUILD_TESTING)
	#Setup CMake to run tests
	

	#Prep for compiling against boost
	find_package(Boost COMPONENTS unit_test_framework REQUIRED)
	include_directories (${Boost_INCLUDE_DIRS})

	add_executable(nups_test test/nups_test.cpp) 
	target_link_libraries(nups_test ${Boost_LIBRARIES} nups)
	add_test(NUPS nups_test)

ENDIF(BUILD_TESTING)







