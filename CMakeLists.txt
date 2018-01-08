#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(Xe1tTools)


#set(CMAKE_INSTALL_PREFIX $ENV{HOME})
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------
#Set the necessary path variables and find the ROOT package

# Tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )

find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore MathMore SPlot Thread Gui)
#find_package(ROOT 5.34 REQUIRED)


set(ROOT_CONFIG_DEBUG ON)

include(${ROOT_USE_FILE})


if(ROOT_CONFIG_DEBUG)
	message(STATUS "CMAKE_PREFIX_PATH: " ${CMAKE_PREFIX_PATH})
	message(STATUS "ROOT_LIBRARY_DIR: " ${ROOT_LIBRARY_DIR})
	foreach(lib ${ROOT_LIBRARIES})
		message(STATUS "Using root library: " ${lib})
	endforeach(lib)
	message(STATUS "ROOT_CXX_FLAGS: " ${ROOT_CXX_FLAGS})
	message(STATUS "ROOT_USE_FILE: " ${ROOT_USE_FILE})
	message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})
endif(ROOT_CONFIG_DEBUG)



SET(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})

#SET(CMAKE_EXE_LINKER_FLAGS "-m64 -L$ENV{ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lGui")
SET(CMAKE_EXE_LINKER_FLAGS "-m64 -pthread -lm -ldl -rdynamic")

include_directories($ENV{HOME}/include)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

link_directories(${ROOT_LIBRARY_DIR})
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
link_directories(${PROJECT_SOURCE_DIR}/RootLib)

#This list is needed for making root dictionaries for the Pax classes that store data in rootfile
#set(pax_header
#	${PROJECT_SOURCE_DIR}/include/PaxTree.hh
#)


set(sources
	${PROJECT_SOURCE_DIR}/src/PaxTreeWrapper.cc
)


#Make a library for the PaxEvent
#add_library(PaxTree SHARED G__PaxTree.cxx)
#add_library(PaxTree SHARED ${PROJECT_SOURCE_DIR}/include/PaxTree.hh)

#ROOT_GENERATE_DICTIONARY(G__PaxTree ${pax_header} MODULE PaxTree LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)





add_executable(MakeMiniTree ${PROJECT_SOURCE_DIR}/prog/MakeMiniTree.cpp ${sources})
target_link_libraries(MakeMiniTree ${ROOT_LIBRARIES} libPaxClasses.so)
#target_link_libraries(MakeMiniTree ${ROOT_LIBRARIES} PaxTree)

