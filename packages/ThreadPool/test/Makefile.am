#@HEADER
# ************************************************************************
# 
#                          ThreadPool Package 
#                 Copyright (2008) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Carter Edwards (hcedwar@sandia.gov) 
# 
# ************************************************************************
#@HEADER

SUBDIRS =

# The following line helps the test harness recover from build errors.

all-local:

include $(top_builddir)/Makefile.export.threadpool

EXEEXT = .exe

noinst_PROGRAMS = test_tpi test_tpi_cpp test_sum 

test_tpi_SOURCES      = test_main.c test_tpi_unit.c test_c_dnax.c test_c_tpi.c test_pthreads.c
test_tpi_DEPENDENCIES = $(top_builddir)/src/libtpi.a
test_tpi_CFLAGS     = $(THREADPOOL_INCLUDES)
test_tpi_LDADD        = $(THREADPOOL_LIBS)

test_tpi_cpp_SOURCES      = test_tpi.cpp
test_tpi_cpp_DEPENDENCIES = $(top_builddir)/src/libtpi.a
test_tpi_cpp_CXXFLAGS     = $(THREADPOOL_INCLUDES)
test_tpi_cpp_LDADD        = $(THREADPOOL_LIBS)

test_sum_SOURCES      = test_mpi_sum.c
test_sum_DEPENDENCIES = $(top_builddir)/src/libtpi.a
test_sum_CFLAGS     = $(THREADPOOL_INCLUDES)
test_sum_LDADD        = $(THREADPOOL_LIBS)

