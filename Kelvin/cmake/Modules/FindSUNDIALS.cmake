#------------------------------------------------------------------------------
# Copyright  2016-, UT-Battelle LLC
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
#
# Find the SUNDIALS library of hybrid CPU+GPU solvers.
#
# Usage:
#   find_package(SUNDIALS [REQUIRED] [QUIET] )
#
# The following variables will be set if the library is found:
#   SUNDIALS_FOUND          ... set to true if module(s) exist
#   SUNDIALS_LIBRARIES      ... list of SUNDIALS libraries (w/o the '-l')
#   SUNDIALS_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
#   SUNDIALS_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
#
# The following variables can be used to influence how the function searches 
# for SUNDIALS:
#   SUNDIALS_ROOT                ... if set, the libraries are exclusively searched
#                                 under this path

# Configure the basic set of libraries
set(SUNDIALS_LIBRARIES_LIST sundials_cvode sundials_nvecserial)

# If Spack is available and SUNDIALS_ROOT specified, see if Spack has an
# installation of Sundials hanging around.
if(NOT SUNDIALS_ROOT AND SPACK_ROOT)
  # Run the Spack search command and deposit the result directly into the
  # SUNDIALS_ROOT variable since it is empty in this branch.
  execute_process(COMMAND "${SPACK_ROOT}/bin/spack" location --install-dir sundials
         OUTPUT_VARIABLE SUNDIALS_ROOT)
  string(STRIP ${SUNDIALS_ROOT} SUNDIALS_ROOT)
endif(NOT SUNDIALS_ROOT AND SPACK_ROOT)

# Debug line. Please leave.
# message(STATUS "SPACK = ${SPACK_ROOT}, SUNDIALS=${SUNDIALS_ROOT}")

# Assign the variables
if(SUNDIALS_ROOT)
  message(STATUS "SUNDIALS found at ${SUNDIALS_ROOT}")
  # Set base paths - this is probably no robust enough if using a system
  # installation on Fedora because it would use lib64 versus lib, although
  # this will work just fine is Spack was used to install SUNDIALS.
  set(SUNDIALS_INCLUDE_DIRS "${SUNDIALS_ROOT}/include")
  set(SUNDIALS_LIBRARY_DIRS "${SUNDIALS_ROOT}/lib")
  # Find each library in the list by full path and add them to the master list.
  foreach(LIB ${SUNDIALS_LIBRARIES_LIST})
       find_library(SUNDIALS_${LIB} NAMES ${LIB} HINTS ${SUNDIALS_LIBRARY_DIRS})
       message(STATUS "Found SUNDIALS library ${LIB} at ${SUNDIALS_${LIB}}.")
       set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARIES} ${SUNDIALS_${LIB}})
  endforeach()
  # Print some information for debugging, etc.
  message(STATUS "Found SUNDIALS libraries: ${SUNDIALS_LIBRARIES}")
  message(STATUS "SUNDIALS Library Directory = ${SUNDIALS_LIBRARY_DIRS}")
  message(STATUS "SUNDIALS Include Directory = ${SUNDIALS_INCLUDE_DIRS}")
  set (SUNDIALS_FOUND TRUE)
elseif(NOT SUNDIALS_ROOT)
  message(FATAL_ERROR "SUNDIALS requested but not found!")
endif(SUNDIALS_ROOT)