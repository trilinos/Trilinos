@HEADER
c ***********************************************************************
c 
c            Trilinos: An Object-Oriented Solver Framework
c                 Copyright (2001) Sandia Corporation
c 
c Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
c license for use of this work by or on behalf of the U.S. Government.
c 
c This library is free software; you can redistribute it and/or modify
c it under the terms of the GNU Lesser General Public License as
c published by the Free Software Foundation; either version 2.1 of the
c License, or (at your option) any later version.
c  
c This library is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c Lesser General Public License for more details.
c  
c You should have received a copy of the GNU Lesser General Public
c License along with this library; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
c USA
c Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
c 
c ***********************************************************************
c @HEADER


c      Dimensions for Fortran test code
c
       INTEGER NZMAX
       INTEGER NUK
       PARAMETER (NUK = 4000,NZMAX = NUK*NUK+2*NUK)
