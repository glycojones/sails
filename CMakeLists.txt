## Project settings
#
cmake_minimum_required(VERSION 3.1)
project(sails C CXX)
option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SAILS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/source/sails)
set(PRIVATEER_DIR ${CMAKE_SOURCE_DIR}/dependencies/privateer)
set(PRIVATEER_SOURCE_DIR ${PRIVATEER_DIR}/src/privateer)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(GSL REQUIRED)

## Find Sails' dependencies that have been compiled locally. 
#
find_library(MMDB2DEP NAMES mmdb2
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CCP4CDEP NAMES ccp4c
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CCP4SRSDEP NAMES ccp4srs
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERCOREDEP NAMES clipper-core
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERMMDBDEP NAMES clipper-mmdb
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERMINIMOLDEP NAMES clipper-minimol
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERCONTRIBDEP NAMES clipper-contrib
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERCCP4DEP NAMES clipper-ccp4
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERCIFDEP NAMES clipper-cif
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(CLIPPERCIFDEP NAMES clipper-cif
            HINTS ${PRIVATEER_DIR}/dependencies/lib)
find_library(PRIVATEERDEP NAMES privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12.1-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12.2-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12.3-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12.4-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12.5-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-12.6-x86_64-3.7/privateer
            HINTS ${PRIVATEER_DIR}/build/lib.macosx-13-x86_64-3.7/privateer)

#Add gemmi that needs to be integrated into sails' shared library
add_subdirectory(${PRIVATEER_DIR}/dependencies/pybind11/)



#Find Sails's main.cpp file
set(SAILS_SOURCE ${SAILS_SOURCE_DIR}/sails.cpp)

#Location of include files to be linked
include_directories(
                    ${PRIVATEER_DIR}/dependencies/include
                    ${PRIVATEER_DIR}/dependencies/gemmi/include
                    ${PRIVATEER_DIR}/dependencies/pybind11/include
                    ${PRIVATEER_DIR}/dependencies/json/include
                    ${PRIVATEER_DIR}/src
                    ${PRIVATEER_DIR}/src/privateer/cpp
                    ${PYTHON_INCLUDE_DIRS}
                    ${SAILS_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/dependencies/include
                    ${CMAKE_SOURCE_DIR}/dependencies/gemmi/include
                    ${PYBIND11_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${PRIVATEER_SOURCE_DIR}/cpp
                    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11)

# ## Create libsails
add_library(sails_lib
            ${SAILS_SOURCE_DIR}/sails-find.cpp
            ${SAILS_SOURCE_DIR}/sails-lib.cpp)


target_link_libraries ( sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY}
                        ${PRIVATEERDEP})


## Create Sails executables
#
add_executable(sails_exec
                ${SAILS_SOURCE})

add_executable(sails_combine_exec
                ${SAILS_SOURCE_DIR}/tools/sails-backend_combine.cpp)

add_executable(sails_fingerprint_exec
                ${SAILS_SOURCE_DIR}/tools/sails-backend_fingerprint.cpp)

add_executable(sails_minimax_exec
                ${SAILS_SOURCE_DIR}/tools/sails-backend_minimax.cpp)             

target_link_libraries ( sails_exec
                        sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY})

target_link_libraries ( sails_combine_exec
                        sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY} )

target_link_libraries ( sails_fingerprint_exec 
                        sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY}
                        GSL::gsl )

target_link_libraries ( sails_minimax_exec     
                        sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY} )

set_target_properties ( sails_lib    PROPERTIES OUTPUT_NAME sails )
set_target_properties ( sails_exec   PROPERTIES OUTPUT_NAME sails )

target_link_libraries ( sails_exec 
                        sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY})

pybind11_add_module(sails_module  ${SAILS_SOURCE_DIR}/sails-pybind11.cpp )
set_target_properties ( sails_module PROPERTIES OUTPUT_NAME sails_core )
target_link_libraries ( sails_module PRIVATE 
                        sails_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY})
                        

install ( TARGETS sails_lib
            LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS sails_exec DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS sails_module DESTINATION ${PROJECT_SOURCE_DIR} )

if (DEVKIT)
  install ( TARGETS sails_combine_exec     DESTINATION bin )
  install ( TARGETS sails_fingerprint_exec DESTINATION bin )
  install ( TARGETS sails_minimax_exec     DESTINATION bin )
endif()
