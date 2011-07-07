#@HEADER
# ************************************************************************
# 
#               Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms
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

#EXTRA_DIST = Makefile.external.stub.in

install-data-local:
	$(mkinstalldirs) $(exec_prefix)/example/thyra ; \
	$(INSTALL_DATA) $(srcdir)/../operator_vector/sillyCgSolve.hpp $(exec_prefix)/example/thyra/. ; \
	$(INSTALL_DATA) $(srcdir)/../operator_vector/sillierCgSolve.hpp $(exec_prefix)/example/thyra/. ; \
	$(INSTALL_DATA) $(srcdir)/../operator_vector/silliestCgSolve.hpp $(exec_prefix)/example/thyra/. ; \
	$(INSTALL_DATA) $(srcdir)/../operator_vector/sillyCgSolve_serial.cpp $(exec_prefix)/example/thyra/. ; \
	$(INSTALL_DATA) $(srcdir)/../operator_vector/ExampleTridiagSerialLinearOp.hpp $(exec_prefix)/example/thyra/. ; \
	$(INSTALL_DATA) Makefile.external.stub $(exec_prefix)/example/thyra/Makefile

sample-makefile:
	cp Makefile.external.stub $(srcdir)/Makefile.external.sample
