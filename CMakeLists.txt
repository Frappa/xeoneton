#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(Xe1tTools)

set(ROOT_CONFIG_DEBUG OFF)
# Tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )

#set(CMAKE_INSTALL_PREFIX $ENV{HOME})
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------
#Set the necessary path variables and find the ROOT package
find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore MathMore RooFitCore RooFit RooStats SPlot Thread Gui)
#find_package(ROOT REQUIRED)


include(${ROOT_USE_FILE})


if(ROOT_CONFIG_DEBUG)
	message(STATUS "ROOT_LIBRARY_DIR: " ${ROOT_LIBRARY_DIR})
	foreach(lib ${ROOT_LIBRARIES})
		message(STATUS "Using root library: " ${lib})
	endforeach(lib)
	message(STATUS "ROOT_CXX_FLAGS: " ${ROOT_CXX_FLAGS})
	message(STATUS "ROOT_USE_FILE: " ${ROOT_USE_FILE})
	message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})
endif(ROOT_CONFIG_DEBUG)



SET(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})

include_directories($ENV{HOME}/include)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

link_directories(${ROOT_LIBRARY_DIR} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


#This list is needed for making root dictionaries for the Pax classes that store data in rootfile
set(pax_header
    ${PROJECT_SOURCE_DIR}/include/PaxTree.hh
)


#Make a library for the PaxEvent
#add_library(PaxTree SHARED G__PaxTree.cxx)

#ROOT_GENERATE_DICTIONARY(G__PaxTree ${pax_header} MODULE PaxTree LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)



add_executable(MakeMiniTree ${PROJECT_SOURCE_DIR}/prog/MakeMiniTree.cpp)
#target_link_libraries(MakeMiniTree ${ROOT_LIBRARIES})
target_link_libraries(MakeMiniTree ${ROOT_LIBRARIES} PaxTree)

