#@HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
# 
# ************************************************************************
#@HEADER
# Scripts in this directory are automatically detected by the test harness
# if set up to build properly below.  Scripts are executed by the test harness.

#np# This Makefile.am basically copies the "autotoolsTestAll" script from the 
#np# source tree to the build tree and gives it a .exe extension.  To set up
#np# multiple scripts, list all scripts in the "PROGRAMS" primary and create
#np# a (script_name)_SOURCES primary and make rule for each script.  Don't 
#np# forget to get rid of the "autotoolsTestAll" example.  This  directory 
#np# is for serial scripts only.  Keep all rules and primaries within 
#np# the "else" portion of the "if HAVE_MPI" statement.

EXEEXT = .exe
# Nothing should be done for an mpi build
if HAVE_MPI
else

noinst_PROGRAMS = autotoolsTestAll


autotoolsTestAll_SOURCES = $(srcdir)/autotoolsTestAll

autotoolsTestAll$(EXEEXT): $(srcdir)/autotoolsTestAll
	cp $(srcdir)/autotoolsTestAll autotoolsTestAll$(EXEEXT) 

endif
