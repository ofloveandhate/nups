#useful defines:

# -DCMAKE_GENERATOR_PLATFORM=x64 or =x86



cmake_minimum_required (VERSION 3.4)
project (nups)

# The version number.
set (nups_VERSION_MAJOR 0)
set (nups_VERSION_MINOR 1)

set(CMAKE_CXX_STANDARD 11)


if (WIN32)
	# disable autolinking in boost
	add_definitions( -DBOOST_ALL_DYN_LINK )
	add_definitions( -DBOOST_ALL_NO_LIB )
	set(Boost_USE_MULTITHREADED 	ON)
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)
endif()	


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


include(GenerateExportHeader)


set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#include_directories("${PROJECT_BINARY_DIR}")

add_library(nups SHARED ${NUPS_HEADERS} ${NUPS_SOURCES})

generate_export_header(nups           
    BASE_NAME nups
    EXPORT_MACRO_NAME NUPS_EXPORTS
    EXPORT_FILE_NAME nups_exports.h
    STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)

target_include_directories(nups PUBLIC ${CMAKE_BINARY_DIR})

install(
		TARGETS nups 
		ARCHIVE DESTINATION "lib"
		LIBRARY DESTINATION "lib"
		COMPONENT library)

install(
		FILES ${NUPS_HEADERS} 
		DESTINATION "include/nups")


enable_testing()
#Setup CMake to run tests



#Prep for compiling against boost
find_package(Boost REQUIRED
			COMPONENTS 
				unit_test_framework )

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ${CMAKE_BINARY_DIR}/lib)


set(NUPS_TEST_SOURCES
	test/nups_test.cpp)

add_executable(nups_test ${NUPS_TEST_SOURCES}) 
target_link_libraries(nups_test ${Boost_LIBRARIES} nups)

add_test(NAME nups_test 
		COMMAND ${CMAKE_BINARY_DIR}/bin/nups_test)









