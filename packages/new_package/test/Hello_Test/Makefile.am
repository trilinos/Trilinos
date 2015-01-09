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

#np# Include the Makefile.export file for this package to provide
#np# access to the _LIBS and _INCLUDES macros
include $(top_builddir)/Makefile.export.new_package

#np# List the names of all exe's in this directory here
bin_PROGRAMS = Hello_Test.exe

#np# One SOURCES primary is needed for EACH exe listed above.  (Notice for all
#np# of the following primaries that the pattern is
#np# "(Exe_name)_(PRIMARY_NAME) = ").
#np# List all source files needed to compile the exe.
Hello_Test_exe_SOURCES = $(srcdir)/Hello_Test.cpp

#np# One DEPENDENCIES primary is used for each exe listed above.  List any
#np# libraries that the exe depends on.  (This primary is not absolutely
#np# essential.)
Hello_Test_exe_DEPENDENCIES = $(top_builddir)/src/libnew_package.a $(top_builddir)/../epetra/src/libepetra.a

#np# One LDADD primary can be defined for each exe listed above.  Include
#np# anything that should be passed to ld.  This is typically taken care 
#np# of by the package's _LIBS value from Makefile.export.<package>
Hello_Test_exe_LDADD = $(NEW_PACKAGE_LIBS)

#np# One CPPFLAGS primary can be defined for each exe listed above.  Include
#np# any -I's necessary for building the exe.  This is typically taken care
#np# of by the package's _INCLUDES value from Makefile.export.<package>
AM_CPPFLAGS = $(NEW_PACKAGE_INCLUDES)

