# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
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

all:
	@echo "Using TRILINOS_HOME=$(TRILINOS_HOME)"
	@echo "Using TRILINOS_ARCH=$(TRILINOS_ARCH)"
	@echo "Using TRILINOS_LIB=$(TRILINOS_LIB)"
	(cd anasazi; make; cd -);
	(cd amesos; make; cd -);
	(cd aztecoo; make; cd -);
	(cd epetra; make; cd -);
	(cd ifpack; make; cd -);
	(cd ml; make; cd -);
	(cd nox; make; cd -);
	(cd triutils; make; cd -);
	@echo "NOTE: Amesos examples have not been complied";

clean:
	(cd amesos; make clean; cd -);
	(cd anasazi; make clean; cd -);
	(cd aztecoo; make clean; cd -);
	(cd epetra; make clean; cd -);
	(cd ifpack; make clean; cd -);
	(cd ml; make clean; cd -);
	(cd nox; make clean; cd -);
	(cd triutils; make clean; cd -);
	(cd triutils; make clean; cd -);
