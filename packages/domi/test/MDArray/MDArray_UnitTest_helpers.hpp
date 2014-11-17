/*
// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Domi_MDArray.hpp"
#include "Domi_MDArrayRCP.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace MDArrayUnitTestHelpers {


extern int nrows;
extern int ncols;
extern int nlevs;


template< class T >
Domi::MDArray< T > generateMDArray(const int nrows_in,
                                   const int ncols_in,
                                   Domi::Layout layout =
                                     Domi::DEFAULT_ORDER)
{
  using Teuchos::tuple;
  using Teuchos::as;

  typedef typename Domi::dim_type dim_type;

  Domi::MDArray< T > a(tuple< dim_type >(nrows_in, ncols_in),
                       layout);
  T value = 0;
  for( int j = 0; j < ncols_in; ++j )
    for( int i = 0; i < nrows_in; ++i )
      a(i,j) = value++; // tests non-const operator()(i,...)
  return a;
}


template< class T >
Domi::MDArray< T > generateMDArray(const int nrows_in,
                                   const int ncols_in,
                                   const int nlevs_in,
                                   Domi::Layout layout =
                                     Domi::DEFAULT_ORDER)
{
  using Teuchos::tuple;
  using Teuchos::as;

  typedef typename Domi::dim_type dim_type;

  Domi::MDArray< T > a(tuple< dim_type >(nrows_in, ncols_in, nlevs_in),
                       layout);
  T value = 0;
  for( int i = 0; i < nrows_in; ++i )
    for( int j = 0; j < ncols_in; ++j )
      for ( int k = 0; k < nlevs_in; ++k)
	a(i,j,k) = value++; // tests non-const operator()(i,...)

  return a;
}


template< class T >
Domi::MDArrayRCP< T > generateMDArrayRCP(const int nrows_in,
                                         const int ncols_in,
                                         Domi::Layout layout =
                                           Domi::DEFAULT_ORDER)
{
  using Teuchos::tuple;
  using Teuchos::as;

  typedef typename Domi::dim_type dim_type;

  Domi::MDArrayRCP< T > a(tuple< dim_type >(nrows_in, ncols_in),
                          layout);
  T value = 0;
  for( int j = 0; j < ncols_in; ++j )
    for( int i = 0; i < nrows_in; ++i )
      a(i,j) = value++; // tests non-const operator()(i,...)
  return a;
}


template< class T >
Domi::MDArrayRCP< T > generateMDArrayRCP(const int nrows_in,
                                         const int ncols_in,
                                         const int nlevs_in,
                                         Domi::Layout layout =
                                           Domi::DEFAULT_ORDER)
{
  using Teuchos::tuple;
  using Teuchos::as;

  typedef typename Domi::dim_type dim_type;

  Domi::MDArrayRCP< T > a(tuple< dim_type >(nrows_in, ncols_in, nlevs_in),
                          layout);
  T value = 0;
  for( int i = 0; i < nrows_in; ++i )
    for( int j = 0; j < ncols_in; ++j )
      for ( int k = 0; k < nlevs_in; ++k)
	a(i,j,k) = value++; // tests non-const operator()(i,...)

  return a;
}


} // namespace ArrayUnitTestHelpers
