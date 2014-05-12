#@HEADER
# ************************************************************************
# 
#                               FEI
#                 Copyright (2007) Sandia Corporation
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
# Questions? Contact Alan B. Williams (william@sandia.gov) 
# 
# ************************************************************************
#@HEADER

DEFS = 
CXXLD=$(CXX)

EXEEXT = .exe

noinst_PROGRAMS = copyinputfiles

copyinputfiles_SOURCES = \
  $(top_srcdir)/test/multifield/mf.i

copyinputfiles$(EXEEXT): $(copyinputfiles_SOURCES)
	@if [ "$(top_builddir)" != "$(top_srcdir)" ] ; then \
	  echo " "; \
	  echo "copying 'src-dir'/test/multifield/*.* to 'build-dir'/test/multifield";\
	  echo " "; \
	  cp $(top_srcdir)/test/multifield/*.* $(top_builddir)/test/multifield; \
	fi

run-tests:copyinputfiles$(EXEEXT)
	$(RUN_COMMAND) 1 ../fei_test.exe -i mf.i

