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
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories("${PROJECT_BINARY_DIR}")

#IF (BUILD_SHARED_LIBS)
#	add_library(nups SHARED ${NUPS_HEADERS} ${NUPS_SOURCES})
#ELSE ()
#	
#ENDIF(BUILD_SHARED_LIBS)

add_library(nups STATIC ${NUPS_HEADERS} ${NUPS_SOURCES})

install(
		TARGETS nups 
		ARCHIVE DESTINATION "lib"
		LIBRARY DESTINATION "lib"
		COMPONENT library)

install(
		FILES ${NUPS_HEADERS} 
		DESTINATION "include/nups")



IF (BUILD_TESTING)
	enable_testing()
	#Setup CMake to run tests
	
	#set(Boost_USE_STATIC_LIBS       ON)
	#
	#add_definitions(-DBOOST_ALL_NO_LIB)  # tell the Boost config system not to automatically select which libraries to link against.
	
	#set(Boost_NO_SYSTEM_PATHS 		true)
	#set(Boost_USE_STATIC_LIBS 		ON)# CACHE BOOL "use static libraries from Boost")
	set(Boost_USE_MULTITHREADED 	ON)
	
	#set(Boost_USE_STATIC_LIBS        OFF)
	#set(Boost_USE_STATIC_RUNTIME     OFF)
	add_definitions( -DBOOST_ALL_NO_LIB )
	add_definitions( -DBOOST_ALL_DYN_LINK )
	
	if (WIN32)
			# disable autolinking in boost
		

		# force all boost libraries to dynamic link (we already disabled
		# autolinking, so I don't know why we need this, but we do!)
		
	endif()	
	
	
	#Prep for compiling against boost
	find_package(Boost REQUIRED
				COMPONENTS 
					unit_test_framework )
	
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
	LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
	
	
	
	
	add_executable(nups_test test/nups_test.cpp) 
	target_link_libraries(nups_test ${Boost_LIBRARIES} nups)
	add_test(NUPS nups_test)

ENDIF(BUILD_TESTING)







