/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifdef CBLAS 
#define XCOPY ccopy
#define XSCAL cscal
#define XAXPY caxpy
#define IXAMAX icamax
#define XASUM scasum 
#define XDOT cdotu

#else

#define XCOPY(len,v1,s1,v2,s2) ccopy_(&len,v1,&s1,v2,&s2)
#define XSCAL(len,scal,v1,s1) cscal_(&len,&scal,v1,&s1)
#define XAXPY(len,scal,v1,s1,v2,s2) caxpy_(&len,&scal,v1,&s1,v2,&s2)
#define IXAMAX(len,v1,s1) icamax_(&len,v1,&s1)
#define XASUM(v3,len,v1,s1) scasum_(v3,&len,v1,&s1)
#define XDOT(v3,len,v1,s1,v2,s2) cdotu_(v3,&len,v1,&s1,v2,&s2)


#endif
