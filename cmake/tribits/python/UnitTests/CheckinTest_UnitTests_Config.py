# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

configuration = {

    'defaults': {
        '--send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
        },

    'cmake': {
        
        'common': [
            '-DTPL_ENABLE_Pthread:BOOL=OFF',
            '-DTPL_ENABLE_BinUtils:BOOL=OFF',
            ],
        'default-builds': [    
          # This is a list of tuples so we can preserve order.
          ('MPI_DEBUG', [
            '-DTPL_ENABLE_MPI:BOOL=ON',
            '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
            '-DTrilinos_ENABLE_DEBUG:BOOL=ON',
            '-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON',
            '-DTrilinos_ENABLE_DEBUG_SYMBOLS:BOOL=ON',
            '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON',
            '-DTeuchos_ENABLE_DEFAULT_STACKTRACE:BOOL=OFF',
            ]),
          # Order matters here for the expected test output.
          ('SERIAL_RELEASE', [
            '-DTPL_ENABLE_MPI:BOOL=OFF',
            '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
            '-DTrilinos_ENABLE_DEBUG:BOOL=OFF',
            '-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF',
            '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF',
            ]),
          ]
    },

}

