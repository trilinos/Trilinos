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

if HAVE_AMESOS_SLUS
SLUS_TEST = TestAmesosSerial
else
SLUS_TEST = 
endif

if HAVE_AMESOS_UMFPACK
UMFPACK_TEST = TestUmfpackSerial
else
UMFPACK_TEST = 
endif

if HAVE_AMESOS_KLU
KLU_TEST = TestKluSerial
else
KLU_TEST = 
endif

noinst_PROGRAMS = $(SLUS_TEST) $(UMFPACK_TEST) $(KLU_TEST)

if HAVE_AMESOS_SLUS
TestAmesosSerial_SOURCES = $(srcdir)/TestAmesosSerial

TestAmesosSerial$(EXEEXT): $(srcdir)/TestAmesosSerial
	cp $(srcdir)/TestAmesosSerial TestAmesosSerial$(EXEEXT) 
endif

if HAVE_AMESOS_UMFPACK
TestUmfpackSerial_SOURCES = $(srcdir)/TestUmfpackSerial

TestUmfpackSerial$(EXEEXT): $(srcdir)/TestUmfpackSerial
	cp $(srcdir)/TestUmfpackSerial TestUmfpackSerial$(EXEEXT) 
endif

if HAVE_AMESOS_KLU
TestKluSerial_SOURCES = $(srcdir)/TestKluSerial

TestKluSerial$(EXEEXT): $(srcdir)/TestKluSerial
	cp $(srcdir)/TestKluSerial TestKluSerial$(EXEEXT) 
endif


