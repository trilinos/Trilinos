#@HEADER
# ************************************************************************
# 
#                     New_Package Example Package
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER
#Builds a simple tester program.  Note that executables must end with a .exe
#to allow certain proprietary operating systems to use them.

#This test should only be built if --enable-newp_swahili is specified

if HAVE_NEWP_SWAHILI

include $(top_builddir)/Makefile.export.new_package

bin_PROGRAMS = Jambo_Test.exe

Jambo_Test_exe_SOURCES = $(srcdir)/Jambo_Test.cpp

Jambo_Test_exe_DEPENDENCIES = $(top_builddir)/src/libnew_package.a $(top_builddir)/../epetra/src/libepetra.a

Jambo_Test_exe_LDADD = $(NEW_PACKAGE_LIBS)

AM_CPPFLAGS = $(NEW_PACKAGE_INCLUDES)

else
endif
