## Project settings
#
cmake_minimum_required(VERSION 3.1)
project(sails C CXX)
option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

if (DEVKIT)
  message(STATUS "Building with local devkit options")
endif()

## Required CCP4 CCP4_LIBRARIES plus GSL
#
find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-mmdb clipper-minimol
clipper-contrib clipper-ccp4 clipper-cif )

find_package(GSL REQUIRED)

## Required Python stuff - link against system Python3 for testing on YSBL machines
#
if (DEVKIT)
    if(APPLE)
        find_program(PYTHON NAMES $ENV{PYTHON} python3.7
                            HINTS ${CMAKE_INSTALL_PREFIX}/bin)
    else()
        find_program(PYTHON NAMES $ENV{PYTHON} python3.6m python3.6 python3
                            HINTS ${CMAKE_INSTALL_PREFIX}/bin
                            HINTS /usr/local/bin
                            HINTS /usr/bin )
    endif()
else()
    # Link against CCP4-distributed Python
    if(APPLE)
        find_program(PYTHON NAMES $ENV{PYTHON} ccp4-python
                            HINTS ${CMAKE_INSTALL_PREFIX}/bin )
    else()
        find_program(PYTHON NAMES $ENV{PYTHON} ccp4-python python.exe
                            HINTS ${CMAKE_INSTALL_PREFIX}/bin
                            PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
    endif()
endif()

## Configure Python
#
execute_process(COMMAND ${PYTHON} -c "import sys; sys.stdout.write(sys.prefix)"
                OUTPUT_VARIABLE _PY_PREFIX
                RESULT_VARIABLE failure)
if (failure)
    message(FATAL_ERROR "Failed to run ${PYTHON}: ${failure}")
else()
    message(STATUS "Python interpreter: ${PYTHON}")
    message(STATUS "Python sys.prefix: ${_PY_PREFIX}")
endif()
# Get Python libraries
find_library(PYTHON_LIBRARY NAMES python3.7 python3.6m python
             HINTS ${_PY_PREFIX}/lib ${_PY_PREFIX}/lib64 ${_PY_PREFIX}/libs
             NO_DEFAULT_PATH)
if (NOT PYTHON_LIBRARY)
    message(FATAL_ERROR "Python library not found.")
else ()
    message(STATUS "Python library: ${PYTHON_LIBRARY}")
    set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})
endif()
# Get Python includes
if(DEVKIT)
  if (APPLE)
    find_path(PYTHON_INCLUDE_DIR NAMES Python.h
              HINTS /usr/local/Cellar/python/3.7.3/Frameworks/Python.framework/Versions/3.7/include/python3.7m
              NO_DEFAULT_PATH )
  else()
    find_path(PYTHON_INCLUDE_DIR NAMES Python.h
              HINTS ${_PY_PREFIX}/include
                    ${_PY_PREFIX}/include/python3.7
                    ${_PY_PREFIX}/include/python3.6m
              NO_DEFAULT_PATH)
  endif()
else()
  find_path(PYTHON_INCLUDE_DIR NAMES Python.h
            HINTS ${_PY_PREFIX}/include
                  ${_PY_PREFIX}/include/python2.7
            NO_DEFAULT_PATH)
endif()

if (NOT PYTHON_INCLUDE_DIR)
    message(FATAL_ERROR "Python header not found.")
else ()
    message(STATUS "Python.h is in: ${PYTHON_INCLUDE_DIR}")
endif()

## Get Pybind11 through cmake
#
find_package(pybind11 CONFIG REQUIRED)
MESSAGE( "-- Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
MESSAGE( "--   pybind11_INCLUDE_DIRS = ${pybind11_INCLUDE_DIRS}" )
MESSAGE( "--   pybind11_LIBRARIES    = ${pybind11_LIBRARIES}" )

include_directories(${CCP4_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${pybind11_INCLUDE_DIRS}
                    "${PYTHON_INCLUDE_DIR}")
## Create libsails
#
add_library(sails_lib
            sails-find.cpp
            sails-lib.cpp)

## Add libprivateer - better mechanism needed
#
add_library ( privateer SHARED IMPORTED )

if (APPLE)
  set_target_properties ( privateer PROPERTIES
                          IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/lib/libprivateer.dylib" )
else()
  set_target_properties ( privateer PROPERTIES
                          IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/lib/libprivateer.so" )
endif()

target_link_libraries ( sails_lib ${CCP4_LIBRARIES} privateer )

## Create Sails executables
#
add_executable(sails_exec
               sails.cpp)

add_executable(sails_combine_exec
               tools/sails-backend_combine.cpp)

add_executable(sails_fingerprint_exec
               tools/sails-backend_fingerprint.cpp)

add_executable(sails_minimax_exec
               tools/sails-backend_minimax.cpp)

target_link_libraries ( sails_exec             sails_lib ${CCP4_LIBRARIES} )
target_link_libraries ( sails_combine_exec     sails_lib ${CCP4_LIBRARIES} )
target_link_libraries ( sails_fingerprint_exec sails_lib ${CCP4_LIBRARIES} GSL::gsl )
target_link_libraries ( sails_minimax_exec     sails_lib ${CCP4_LIBRARIES} )


set_target_properties ( sails_lib    PROPERTIES OUTPUT_NAME sails )
set_target_properties ( sails_exec   PROPERTIES OUTPUT_NAME sails )

## Make hybrid Python module
#
pybind11_add_module(sails_module sails-pybind11.cpp)
set_target_properties ( sails_module PROPERTIES OUTPUT_NAME sails_core )
target_link_libraries ( sails_module PRIVATE sails_lib ${CCP4_LIBRARIES} ${pybind11_LIBRARIES} )

## And finally, install everything in their reet places
#
set     ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/py3 )
message ( STATUS "Installing Python module to: ${PYTHON_MODULE_PATH}")
file (COPY "python/sails" DESTINATION ${PYTHON_MODULE_PATH})
install(FILES sails-lib.h
        DESTINATION include/sails)
install(TARGETS sails_lib
        LIBRARY DESTINATION "${LIB_INSTALL_DIR}")
install ( TARGETS sails_exec DESTINATION bin)
install ( TARGETS sails_module DESTINATION ${PYTHON_MODULE_PATH}/sails )

if (DEVKIT)
  install ( TARGETS sails_combine_exec     DESTINATION bin )
  install ( TARGETS sails_fingerprint_exec DESTINATION bin )
  install ( TARGETS sails_minimax_exec     DESTINATION bin )
endif()
