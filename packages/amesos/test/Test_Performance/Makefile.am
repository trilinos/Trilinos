#@HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

#
#  All of the AmeosShellScript and TestAmesosScript stuff is designed to 
#  copy some files form the source tree to the build tree
#
#  The idea behind the mechanism is to fool autotools into thinking that
#  AmesosShellScript is an executable and that the way that you create
#  that executable is the "cp" command.
#
#  The files are renamed as they are copied so that configure can work 
#  within the source directory (even though we recommend that users 
#  run configure from a separate build directory.
#

noinst_PROGRAMS = PerfTest \
                  Test_UmfpackPerformance

Test_UmfpackPerformance_SOURCES = $(srcdir)/In_Test_UmfpackPerformance.csh

Test_UmfpackPerformance$(EXEEXT): $(srcdir)/In_Test_UmfpackPerformance.csh 
	cp  $(srcdir)/In_Test_UmfpackPerformance.csh Test_UmfpackPerformance.csh

Matrices:  ${tac_with_trilinosdatadir}/amesos/662_bus.rsa.gz 
	gunzip -c ${tac_with_trilinosdatadir}/amesos/662_bus.rsa.gz > 662_bus.rsa

SRC_CXX = 

PerfTest_SOURCES = $(srcdir)/Amesos_Performance.cpp \
	 $(SRC_CXX)

if HAVE_AMESOS_EPETRAEXT
   EPETRA_EXT_INCLUDE =	-I$(top_srcdir)/../epetraext/src/transform -I$(top_builddir)/../epetraext/src/transform -I$(top_srcdir)/../epetraext/src -I$(top_builddir)/../epetraext/src 
else
   EPETRA_EXT_INCLUDE = 
endif 

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src \
	${EPETRA_EXT_INCLUDE} \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src 

if HAVE_AMESOS_EPETRAEXT
   EPETRA_EXT_LIB = 	$(top_builddir)/../epetraext/src/libepetraext.a 
else
   EPETRA_EXT_LIB = 
endif 

PerfTest_DEPENDENCIES = \
	${EPETRA_EXT_LIB} \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

PerfTest_CXXFLAGS = -I$(top_srcdir)/src

#
#  libamesos.a is loaded twice.  The second time allows umfpack to pick up AMD.
#
PerfTest_LDADD = $(top_builddir)/src/libamesos.a \
	${EPETRA_EXT_LIB} \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(AMESOS_LIBS)  $(top_builddir)/src/libamesos.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS) $(FLIBS)

EXTRA_DIST = $(srcdir)/662_bus.rsa 

clean:  clean-am
	rm -f  662_bus_out.rsa 


