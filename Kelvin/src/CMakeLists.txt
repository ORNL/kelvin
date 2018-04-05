#------------------------------------------------------------------------------
# Copyright 2018-, UT-Battelle, LLC
# All rights reserved.
#
# Author Contact: Jay Jay Billings, jayjaybillings <at> gmail <dot> com
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

#Set the package name
SET(PACKAGE_NAME "KELVIN")
#Set the description
SET(PACKAGE_DESCRIPTION "Kelvin Library")
#Set the library name
SET(LIBRARY_NAME "Kelvin")

#Log that this project will be built
MESSAGE(STATUS "----- Detected and building ${PACKAGE_NAME} -----")

#Look for Parsers and MFEM
if (PARSERS_DIR AND MFEM_DIR)

   #Be nice
   message(STATUS "Parsers and MFEM root directories set. Attempting to link.")
   
   #Link to mfem. Note that the value of MFEM_DIR will be overwritten by 
   #find_package!
   find_package(MFEM CONFIG HINTS ${MFEM_DIR}/lib/cmake/mfem)

   # Add the variables to the global property list
   set(${PACKAGE_NAME}_LIBRARY_DIRS ${MFEM_LIBRARY_DIR} CACHE INTERNAL "${PACKAGE_NAME}_LIBRARY_DIRS")
   set(${PACKAGE_NAME}_LIBRARIES ${MFEM_LIBRARY_DIR}/lib${MFEM_LIBRARIES}.a ${LIBRARY_NAME} CACHE INTERNAL "${PACKAGE_NAME}_LIBRARIES")
   set(${PACKAGE_NAME}_INCLUDE_DIRS ${PARSERS_DIR}/include/ ${MFEM_INCLUDE_DIRS} CACHE INTERNAL "${PACKAGE_NAME}_INCLUDE_DIRS")

   #Collect all header filenames in this project 
   #and glob them in HEADERS
   file(GLOB HEADERS *.h)

   #Grab all of the source files
   file(GLOB SRC *.cpp)

   #Add the source code to the library
   add_library(${LIBRARY_NAME} STATIC ${SRC})
   #Link to parsers
   find_library(MFEM_LIBRARY NAMES libmfem.a mfem HINTS ${MFEM_LIBRARY_DIR})
   target_link_libraries(${LIBRARY_NAME} ${MFEM_LIBRARY})
   target_include_directories(${LIBRARY_NAME} PUBLIC ${PARSERS_DIR}/include/ ${MFEM_INCLUDE_DIRS})
    
   #Get the test files
   file(GLOB test_files tests/*Test.cpp)
   #Add the tests
   
      #FIXME! ADD INCLUDE DIRECTORIES then pass to test
   
   set(all_includes_dirs ${CMAKE_CURRENT_SOURCE_DIR} ${PARSERS_DIR}/include ${MFEM_INCLUDE_DIRS})
 #  add_tests("${test_files}" "${all_includes_dirs}" "${MFEM_LIBRARIES}")   

   #Install the Fire header files   
   install(FILES ${HEADERS} DESTINATION include)

else ()

   #Complain
   message(STATUS "Parsers and/or MFEM not found. Please set their \
install directories.")

endif (PARSERS_DIR AND MFEM_DIR)
