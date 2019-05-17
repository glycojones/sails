cmake_minimum_required(VERSION 3.1)
project(sails C CXX)

option(BUILD_SHARED_LIBS "Build as a shared library" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

find_package(CCP4 REQUIRED mmdb2 ccp4c ccp4srs clipper-core clipper-mmdb clipper-minimol
clipper-contrib clipper-ccp4 clipper-cif )

find_package(GSL REQUIRED)

## Find Python stuff
#

if(APPLE)
    find_program(PYTHON NAMES $ENV{PYTHON} python3.7 ccp4-python
                        HINTS /usr/local/bin /usr/bin ${CMAKE_INSTALL_PREFIX}/bin )
else()
    find_program(PYTHON NAMES $ENV{PYTHON} python3.7 python3.6 ccp4-python python.exe
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin
                        PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
endif()


execute_process(COMMAND ${PYTHON} -c "import sys; sys.stdout.write(sys.prefix)"
                OUTPUT_VARIABLE _PY_PREFIX
                RESULT_VARIABLE failure)

if (failure)
    message(FATAL_ERROR "Failed to run ${PYTHON}: ${failure}")
else()
    message(STATUS "Python interpreter: ${PYTHON}")
endif()

message(STATUS "Python sys.prefix: ${_PY_PREFIX}")

if (APPLE)
  find_library(PYTHON_LIBRARY NAMES python3.7 python
               HINTS /usr/local/lib ${_PY_PREFIX}/lib
               NO_DEFAULT_PATH)
else()
  find_library(PYTHON_LIBRARY NAMES python3.7 python3.6m python
               HINTS ${_PY_PREFIX}/lib64 ${_PY_PREFIX}/lib ${_PY_PREFIX}/libs
               NO_DEFAULT_PATH)
endif()

if (NOT PYTHON_LIBRARY)
    message(FATAL_ERROR "Python library not found.")
endif()

message(STATUS "Python library: ${PYTHON_LIBRARY}")

# Rely on Homebrew Python 3.7 until CCP4 bundle it in
# TODO: update upon release

if (APPLE)
find_path(PYTHON_INCLUDE_DIR NAMES Python.h
          HINTS /usr/local/Frameworks/Python.framework/Headers
                ${_PY_PREFIX}/include/python
          NO_DEFAULT_PATH)
else()
  find_path(PYTHON_INCLUDE_DIR NAMES Python.h
            HINTS ${_PY_PREFIX}/include/python3.6m
            NO_DEFAULT_PATH)
endif()

message(STATUS "Python.h is in: ${PYTHON_INCLUDE_DIR}")

if (NOT PYTHON_INCLUDE_DIR)
    message(FATAL_ERROR "Python header not found.")
endif()

set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})


## OpenMP?
#  Fine for now, but we are using C++11 threads on CCP4 7.1

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(sails_module sails-pybind11.cpp)

MESSAGE( "-- Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
MESSAGE( "--   pybind11_INCLUDE_DIRS = ${pybind11_INCLUDE_DIRS}" )
MESSAGE( "--   pybind11_LIBRARIES    = ${pybind11_LIBRARIES}" )

include_directories(${CCP4_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${pybind11_INCLUDE_DIRS}
                    "${PYTHON_INCLUDE_DIR}")

add_library(sails_lib
            sails-find.cpp
            sails-lib.cpp)

## Add libprivateer
#   This will change in the future. Will produce a suitable cmake file for it

add_library ( privateer SHARED IMPORTED )

if (APPLE)
  set_target_properties ( privateer PROPERTIES
                          IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/lib/libprivateer.dylib" )
else()
  set_target_properties ( privateer PROPERTIES
                          IMPORTED_LOCATION "${CMAKE_INSTALL_PREFIX}/lib/libprivateer.so" )
endif()

target_link_libraries ( sails_lib ${CCP4_LIBRARIES} privateer )

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
set_target_properties ( sails_module PROPERTIES OUTPUT_NAME sails )

target_link_libraries ( sails_module PRIVATE sails_lib ${CCP4_LIBRARIES} ${pybind11_LIBRARIES} )

set ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/py2 )

message(STATUS "Python module path: ${PYTHON_MODULE_PATH}")


install(FILES sails-lib.h
        DESTINATION include/sails)
install(TARGETS sails_lib
        LIBRARY DESTINATION "${LIB_INSTALL_DIR}")
install ( TARGETS sails_exec DESTINATION bin)
install ( TARGETS sails_module DESTINATION ${PYTHON_MODULE_PATH} )


if (UNIX)
  install(CODE "
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink sails
    \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sails\")
        ")
endif()
