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
# * Neither the name of Kelvin nor the names of its
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

# Set the package name
SET(PACKAGE_NAME "PMGen")
# Set the description
SET(PACKAGE_DESCRIPTION "Particle Mesh Generator")

# Log that this project will be built
MESSAGE(STATUS "----- Detected and building ${PACKAGE_NAME} -----")

# Add the Modules directory to pick up extra *.cmake files.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Include the testing framework
include(tests)

# Include directories for the main modules. Note that include 
# directories must come before add subdirectory!
include_directories("${CMAKE_SOURCE_DIR}/src")

# Add the executable
add_executable(PMGen PMGen.cpp)
target_link_libraries(PMGen ${Kelvin_LIBRARIES})
target_include_directories(PMGen PUBLIC ${Kelvin_INCLUDE_DIRS})

# Install the executable
install(TARGETS PMGen RUNTIME DESTINATION bin)