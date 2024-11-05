# @HEADER
# *****************************************************************************
#           Trilinos: An Object-Oriented Solver Framework
#
# Copyright 2001-2024 NTESS and the Trilinos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# A) Start out by just finding the headers and libraries
#




tribits_tpl_find_include_dirs_and_libraries( BinUtils
  REQUIRED_HEADERS link.h bfd.h
  MUST_FIND_ALL_HEADERS
  REQUIRED_LIBS_NAMES bfd iberty
  MUST_FIND_ALL_LIBS
  NO_PRINT_ENABLE_SUCCESS_FAIL
  )


#
# B) Now make sure that you can actually link a program
#


include(CheckCXXSourceCompiles)
include(MultilineSet)
include(PrintVar)


function(check_for_binutils_stacktrace  VARNAME)

  set(SOURCE
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
  
  set(CMAKE_REQUIRED_INCLUDES ${TPL_BinUtils_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${TPL_BinUtils_LIBRARIES})
  check_cxx_source_compiles("${SOURCE}" ${VARNAME})
  
endfunction()


if (TPL_ENABLE_BinUtils)

  check_for_binutils_stacktrace(HAS_TPL_BINUNTILS_STACKTRACE)

  if (HAS_TPL_BINUNTILS_STACKTRACE)
    message(STATUS "Extended attempt to enable tentatively enabled TPL 'BinUtils' passed!")
  else()
    message(STATUS "Extended attempt to enable tentatively enabled TPL 'BinUtils' failed!  Setting TPL_ENABLE_BinUtils=OFF")
    set(TPL_ENABLE_BinUtils OFF CACHE STRING "autoset" FORCE)
  endif()

endif()


print_var(TPL_ENABLE_BinUtils)
