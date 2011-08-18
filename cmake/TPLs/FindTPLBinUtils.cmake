# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


#
# A) Start out by just finding the headers and libraries
#


INCLUDE(TPLDeclareLibraries)


TPL_DECLARE_LIBRARIES( BinUtils
  REQUIRED_HEADERS link.h bfd.h
  MUST_FIND_ALL_HEADERS
  REQUIRED_LIBS_NAMES bfd iberty
  MUST_FIND_ALL_LIBS
  NO_PRINT_ENABLE_SUCCESS_FAIL
  )


#
# B) Now make sure that you can actually link a program
#


INCLUDE(CheckCXXSourceCompiles)
INCLUDE(MultilineSet)
INCLUDE(PrintVar)


FUNCTION(CHECK_FOR_BINUTILS_STACKTRACE  VARNAME)

  SET(SOURCE
  "
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdarg>
#include <csignal>
#include <execinfo.h>
#include <cxxabi.h>
#include <link.h>
#include <bfd.h>


int shared_lib_callback(struct dl_phdr_info *info,
  size_t size, void *_data)
{
  (void)info;
  (void)size;
  (void)_data;
  return 0;
}


void process_section(bfd *abfd, asection *section, void *_data)
{}


void call_a_bunch_of_functions()
{

  void **stacktrace_array = 0;
  const size_t stacktrace_size = backtrace(stacktrace_array, 0);

  int status_demangle = 0;
  char *filename_blob = 0;
  char *d = abi::__cxa_demangle(filename_blob, 0, 0, &status_demangle);

  asection *section = 0;
  bfd *abfd = 0;

  const int status_bfd_gsf = (bfd_get_section_flags(abfd, section) & SEC_ALLOC);

  bfd_vma section_vma = bfd_get_section_vma(abfd, section);

  bfd_size_type section_size = bfd_section_size(abfd, section);

  bfd_vma offset;
  asymbol **symbol_table = 0;
  unsigned int line;

  const char *filename=NULL, *function_name=NULL;
  const int line_found = bfd_find_nearest_line(abfd, section, symbol_table,
    offset, &filename, &function_name, &line);

  const int status_bfd_gff = (bfd_get_file_flags(abfd) & HAS_SYMS);
  void **symbol_table_ptr = 0;
  long n_symbols;
  unsigned int symbol_size;
  n_symbols = bfd_read_minisymbols(abfd, false, symbol_table_ptr, &symbol_size);
  
  abfd = bfd_openr(filename_blob, NULL);

  const int status_bfd_cf = bfd_check_format(abfd, bfd_archive);

  char **matching;
  const int status_bfd_cfm = bfd_check_format_matches(abfd, bfd_object, &matching);

  void *data = 0;
  bfd_map_over_sections(abfd, process_section, &data);
  
  bfd_close(abfd);

  struct match_data *match;
  const int status_dl_iphdr = dl_iterate_phdr(shared_lib_callback, &(*match));
  
}

 
int main()
{
  call_a_bunch_of_functions();
  return 0;
}
"
  )
  
  SET(CMAKE_REQUIRED_INCLUDES ${TPL_BinUtils_INCLUDE_DIRS})
  SET(CMAKE_REQUIRED_LIBRARIES ${TPL_BinUtils_LIBRARIES})
  CHECK_CXX_SOURCE_COMPILES("${SOURCE}" ${VARNAME})
  
ENDFUNCTION()


IF (TPL_ENABLE_BinUtils)

  CHECK_FOR_BINUTILS_STACKTRACE(HAS_TPL_BINUNTILS_STACKTRACE)

  IF (HAS_TPL_BINUNTILS_STACKTRACE)
    MESSAGE(STATUS "Extended attempt to enable tentatively enabled TPL 'BinUtils' passed!")
  ELSE()
    MESSAGE(STATUS "Extended attempt to enable tentatively enabled TPL 'BinUtils' failed!  Setting TPL_ENABLE_BinUtils=OFF")
    SET(TPL_ENABLE_BinUtils OFF CACHE STRING "autoset" FORCE)
  ENDIF()

ENDIF()


PRINT_VAR(TPL_ENABLE_BinUtils)
