#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
#                 Copyright (2005) Sandia Corporation
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = testLOCA

testLOCA_SOURCES = $(srcdir)/testLOCA.py

testLOCA$(EXEEXT): $(srcdir)/testLOCA.py
	cp $(srcdir)/testLOCA.py testLOCA$(EXEEXT)
	chmod +x testLOCA$(EXEEXT)

# ------------------------------------------------------------------------
# For the python modules to be available, we first have to do a 
# 'make install' (or more specifically, install-exec) in the toplevel build
# directory.  We do this using the DESTDIR make macro which preprends 
# DESTDIR to the install prefix.
#
# The dummy_file is a trick to prohibit the 'make install-exec' from 
# running again when it reaches this directory (and thereby creating an
# infinite loop). 
# ------------------------------------------------------------------------

LOCAL_INSTALL_DIR=$(abs_builddir)/install

all-local: dummy_file

dummy_file: 
	touch dummy_file
	cd $(top_builddir) && $(MAKE) DESTDIR=$(LOCAL_INSTALL_DIR) install-exec
	rm -f dummy_file

CLEANFILES = testLOCA$(EXEEXT)

DISTCLEANFILES = setpath.pyc

clean-local:
	rm -rf install
