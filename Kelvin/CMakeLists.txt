#------------------------------------------------------------------------------
# Copyright 2018-, UT-Battelle, LLC
# All rights reserved.
#
# Author Contact: Jay Jay Billings, billingsjj <at> ornl <dot> gov
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of fern nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Author(s): Jay Jay Billings
#----------------------------------------------------------------------------*/

# Set a name for this project.
PROJECT ( Kelvin C CXX )

#Log that this project will be built
MESSAGE(STATUS "----- Configuring Main Kelvin Build -----")

# Set the minimum required version of cmake for this project.
cmake_minimum_required(VERSION 2.8)
# Configure some default options
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
message(STATUS "C++ version ${CXX_STANDARD} configured.")

# Add the Modules directory to pick up extra *.cmake files.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Include the testing framework
include(tests)

# Setup flag overrides if Kelvin is building on its own. This works by 
# determining if a parent is available and allowing flag overrides if not.
get_directory_property(hasParent PARENT_DIRECTORY)
if(NOT hasParent STREQUAL "")
  message(STATUS "Parent directory is ${hasParent}")
  list( APPEND CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{CXX_FLAGS}")
  message(STATUS "Parent detected. Disabling Kelvin flag overrides.")
else()
  message(STATUS "No parent directory found. Stand-alone build assumed. Build flag overrides enabled.")
endif()

# Setup Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.cmake ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" 
    VERBATIM)
endif(DOXYGEN_FOUND)

# Include directories for the main modules. Note that include 
# directories must come before add subdirectory!
include_directories("${CMAKE_SOURCE_DIR}/src")

# Add the main source directories
add_subdirectory(src)

# Add the executable for the original code
add_executable(KelvinTH Kelvin.cpp)
target_link_libraries(KelvinTH ${Kelvin_LIBRARIES})
target_include_directories(KelvinTH PUBLIC ${Kelvin_INCLUDE_DIRS})

#Add the executable for the MPM code
#add_executable(KelvinMPM KelvinMPM.cpp)
#target_link_libraries(KelvinMPM ${KELVIN_LIBRARIES})
#target_include_directories(KelvinMPM PUBLIC ${KELVIN_INCLUDE_DIRS})

# Install the executable
#install(TARGETS myExe mySharedLib myStaticLib
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib/static) 
