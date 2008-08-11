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
# Nothing should be done for a serial build


if HAVE_AMESOS_SUPERLU
SUPERLU_TEST = TestSuperluSerial
TestSuperluSerial_SOURCES = $(srcdir)/TestSuperluSerial

TestSuperluSerial$(EXEEXT): $(srcdir)/TestSuperluSerial
	cp $(srcdir)/TestSuperluSerial TestSuperluSerial$(EXEEXT) 
else
SUPERLU_TEST = 
endif

if HAVE_AMESOS_KLU
KLU_TEST = TestKluSerial
TestKluSerial_SOURCES = $(srcdir)/TestKluSerial

TestKluSerial$(EXEEXT): $(srcdir)/TestKluSerial
	cp $(srcdir)/TestKluSerial TestKluSerial$(EXEEXT) 
else
KLU_TEST = 
endif

if HAVE_AMESOS_UMFPACK
UMFPACK_TEST = TestUmfpackSerial
TestUmfpackSerial_SOURCES = $(srcdir)/TestUmfpackSerial

TestUmfpackSerial$(EXEEXT): $(srcdir)/TestUmfpackSerial
	cp $(srcdir)/TestUmfpackSerial TestUmfpackSerial$(EXEEXT) 
else
UMFPACK_TEST = 
endif

#TestAmesos_SOURCES = $(srcdir)/TestAmesos.sh
#TestAmesos$(EXEEXT):
#	cp $(srcdir)/TestAmesos.sh TestAmesos$(EXEEXT) 

# noinst_PROGRAMS = TestAmesos $(UMFPACK_TEST) $(KLU_TEST)  $(SUPERLU_TEST)
noinst_PROGRAMS = $(UMFPACK_TEST) $(KLU_TEST)  $(SUPERLU_TEST)



