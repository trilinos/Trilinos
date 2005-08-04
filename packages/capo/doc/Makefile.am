# @HEADER
# ************************************************************************
# 
#                            Capo Package
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
# Questions? Contact Joseph Simonis (jpsimon@sandia.gov), 
#	Sandia National Laboratories.
# 
# ************************************************************************
# @HEADER

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = Doxyfile \
	$(srcdir)/images/gold-rule-500x11.gif \
	$(srcdir)/images/gold-rule-557x11.gif \
	$(srcdir)/images/red-vertical-bar-background.gif \
	$(srcdir)/images/snl-globe-left.gif \
	$(srcdir)/images/spacer-16x1.gif \
	$(srcdir)/images/spacer-6x1.gif \
	$(srcdir)/includes/footer.html \
	$(srcdir)/includes/header.html \
	$(srcdir)/includes/stylesheet.css

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

dox: Doxyfile
	cp Doxyfile $(top_srcdir)/doc; cd $(top_srcdir)/doc; doxygen


