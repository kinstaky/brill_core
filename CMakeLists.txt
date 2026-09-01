# Require the cmake minimum version.
cmake_minimum_required(VERSION 3.10)

include(FetchContent)


# Set project properties.
project(
	brill
	VERSION 1.0.0
	DESCRIPTION "A project for analysis RIBLL2026 experiment data."
	LANGUAGES CXX C
)


# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
	message(FATAL_ERROR "You cannot build in a source directory (or any directory with"
						"CmakeLists.txt file). Please make a build subdirectory. Feel free"
						"to remove CMakeCache.txt and CMakeFiles.")
endif()

# Set the build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()


# Only do these if this is the main project, and not if it is included
# through add_subdirectory
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

	# Let's ensure -std=c++xx instead of -std=g++xx
	set(CMAKE_CXX_EXTENSIONS OFF)

	# Set the default cxx standard
	set(CMAKE_CXX_STANDARD 17)

	# Set default compile options
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

	# Let's nicely support folders in IDE's
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)

endif()

# find root
find_package(ROOT 6.20 CONFIG REQUIRED COMPONENTS Core RIO Tree Hist Graf Minuit GenVector MathCore)

# find catima
find_package(catima CONFIG QUIET)
if (NOT catima_FOUND)
	message(STATUS "catima not found locally, fetching it...")
	FetchContent_Declare(
		catima
		GIT_REPOSITORY https://github.com/kinstaky/catima.git
		GIT_TAG v1.90
		GIT_SHALLOW TRUE
		GIT_SHALLOW_DEPTH 1
	)
	FetchContent_MakeAvailable(catima)
endif()

if (TARGET catima::catima)
	set(BRILL_CATIMA_TARGET catima::catima)
elseif (TARGET catima)
	set(BRILL_CATIMA_TARGET catima)
else()
	message(FATAL_ERROR "catima target not found after package resolution.")
endif()


# show properties
MESSAGE("Build Properties")
MESSAGE("  c++ compiler.........: " ${CMAKE_CXX_COMPILER})
MESSAGE("  c++ standard.........: " ${CMAKE_CXX_STANDARD})
MESSAGE("  cmake................: " ${CMAKE_VERSION})
MESSAGE("  generator............: " ${CMAKE_GENERATOR})
MESSAGE("  cern ROOT version....: " ${ROOT_VERSION})
MESSAGE("  cern ROOT path.......: " ${CMAKE_FIND_ROOT_PATH})

add_subdirectory(src)
add_subdirectory(bin)

# Examples are only part of standalone brill builds, not consuming projects.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	add_subdirectory(examples)
endif()

add_library(brill::core ALIAS brill_core)
add_library(brill::t0 ALIAS brill_t0_dssd)
