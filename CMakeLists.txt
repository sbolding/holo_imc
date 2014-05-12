#have to tell cmake the least version it can run, this line must be here
cmake_minimum_required(VERSION 2.6) 
#name of project
project (holo)

#Must have this line to run in cygqin
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

#set to include c++ 11
set(CMAKE_CXX_FLAGS "-std=c++11")

#GLOB will make a list of all the files in src/*.cpp and add them
# to a variables called SOURCES
file(GLOB SOURCES src/*.cpp)

#include_directories(include) #if you want to add a director include to look in, but would need a cmake file in there

#includ
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)


#Create executable based on all the source files in SOURCES
#braces just expand the variable to it's value
add_executable(holo_exec ${SOURCES}) 

#message(STATUS ${CMAKE_SOURCE_DIR}) #print out the source directory, example

#for some reason it ignores the runtime destination command, this
#may just be because it is set above or a cygqin issue
install(TARGETS holo_exec RUNTIME DESTINATION build)


