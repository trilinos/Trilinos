/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/


#include <utility>
#include <iostream>

std::pair<int,int> map_a_to_ab( int n , int i , int j )
{
  std::pair<int,int> result ;

  const int nb = ( n + 1 ) / 2 ;
  const int diag = j - i ;

  if ( diag <= - nb || ( 0 <= diag && diag < nb ) ) {
    result.first = i ;
    result.second = ( n + diag ) % n ;
  }
  else {
    result.first = j ;
    result.second = ( n - diag ) % n ;
  }

  return result ;
}

std::pair<int,int> map_ab_to_a( int n , int i , int k )
{
  std::pair<int,int> result ;

  const int nb = ( n + 1 ) / 2 ;

  result.first = i ;
  result.second = ( i + k ) % n ;

  return result ;
}

int main()
{
  const int n = 5 ;

  for ( int j = 0 ; j < n ; ++j ) {
    for ( int i = 0 ; i < n ; ++i ) {
      std::pair<int,int> y = map_a_to_ab( n , j , i );

      std::cout << " (" << y.first << "," << y.second << ")" ;
    }
    std::cout << std::endl ;
  }

  std::cout << std::endl ;

  for ( int j = 0 ; j < n ; ++j ) {
    for ( int i = 0 ; i < n ; ++i ) {
      std::pair<int,int> x = map_a_to_ab( n , i , j );

      std::cout << " (" << x.first << "," << x.second << ")" ;
    }
    std::cout << std::endl ;
  }

  std::cout << std::endl ;

  const int nb = ( n + 1 ) / 2 ;

  for ( int j = 0 ; j < n ; ++j ) {
    for ( int i = 0 ; i < nb ; ++i ) {
      std::pair<int,int> x = map_ab_to_a( n , j , i );
      std::cout << " (" << x.first << "," << x.second << ")" ;
    }
    std::cout << std::endl ;
  }
}

