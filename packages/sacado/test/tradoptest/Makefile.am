## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

EXEEXT = .exe

noinst_PROGRAMS = \
	tradoptest1 \
	tradoptest2 \
	tradoptest3 \
	tradoptest4 \
	tradoptest5 \
	tradoptest6

testall : $(noinst_PROGRAMS)
	./tradoptest1$(EXEEXT)
	./tradoptest2$(EXEEXT)
	./tradoptest3$(EXEEXT)	
	./tradoptest4$(EXEEXT)
	./tradoptest5$(EXEEXT)
	./tradoptest6$(EXEEXT)

tradoptest1_SOURCES = tradoptest.cpp 
tradoptest1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

tradoptest2_SOURCES = tradoptest.cpp 
tradoptest2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

tradoptest3_SOURCES = tradoptest.cpp 
tradoptest3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

tradoptest4_SOURCES = tradoptest.cpp 
tradoptest4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

tradoptest5_SOURCES = tradoptest.cpp 
tradoptest5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

tradoptest6_SOURCES = tradoptest.cpp 
tradoptest6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

tradoptest.cpp : testgen.as tradoptest.c0
	$(SHELL) testgen.as $(srcdir)/tradoptest.c0 > tradoptest.cpp

CLEANFILES = tradoptest.cpp