/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP
#define TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>

#define IFPACK2_INST(CLASSNAME,S,LO,GO) \
  template class CLASSNAME<Tpetra::CrsMatrix<S, LO, GO> >;

#define IFPACK2_INST_GRAPH(CLASSNAME,LO,GO) \
  template class CLASSNAME<Tpetra::CrsGraph<LO, GO> >;

#define IFPACK2_CLASS_CrsMatrix_float_int_int_defaultNode_defaultOps(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,float,int,int)

#define IFPACK2_CLASS_CrsMatrix_float_short_int_defaultNode_defaultOps(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,float,short,int)

#define IFPACK2_CLASS_CrsMatrix_double_int_int_defaultNode_defaultOps(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,double,int,int)

#define IFPACK2_INSTANT_CRSMATRIX_FLOAT_DOUBLE_DEFAULTS(CLASSNAME) \
  IFPACK2_CLASS_CrsMatrix_double_int_int_defaultNode_defaultOps(CLASSNAME)

#define IFPACK2_INSTANT_CRSMATRIX_COMPLEX_DEFAULTS(CLASSNAME) \
  IFPACK2_INST(CLASSNAME,std::complex<double>,int,int) \
  IFPACK2_INST(CLASSNAME,std::complex<float>,int,int)

#endif // TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP

