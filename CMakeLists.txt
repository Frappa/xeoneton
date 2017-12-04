#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(Xe1tTools)

set(ROOT_CONFIG_DEBUG OFF)
# Tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )

set(CMAKE_INSTALL_PREFIX $ENV{HOME})
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

#----------------------------------------------------------------------------
#Set the necessary path variables and find the ROOT package
find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore MathMore Minuit2 RooFitCore RooFit RooStats SPlot Thread Gui)


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
include_directories($ENV{ROOTSYS}/include)
include_directories(${PROJECT_SOURCE_DIR}/include)

link_directories(${ROOT_LIBRARY_DIR})

#This list is needed for making root dictionaries
set(headers
	${PROJECT_SOURCE_DIR}/include/PaxTree.hh
)


set(dict_opt ${dict_opt} "-fPIC")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  #message(STATUS "dir='${dir}'")
  set(dict_opt ${dict_opt} "-I${dir}")
endforeach()
string(REPLACE ";" " " dict_opt "${dict_opt}")


#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)


#message(STATUS "Options for dictionary generation = ${dict_opt}")
ROOT_GENERATE_DICTIONARY(G__PaxStruct ${headers} MODULE PaxStruct LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)


add_executable(MakeMiniTree ${PROJECT_SOURCE_DIR}/prog/MakeMiniTree.cpp G__PaxStruct.cxx);
target_link_libraries(MakeMiniTree ${ROOT_LIBRARIES})

