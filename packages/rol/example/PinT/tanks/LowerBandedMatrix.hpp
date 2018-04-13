// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef LOWER_BANDED_MATRIX_HPP
#define LOWER_BANDED_MATRIX_HPP

#include<algorithm>
#include<vector>
#include<exception>
#include<iostream>
#include<iomanip>

/** \class LowerBandedMatrix
    \brief Implements GAXPY-like operations with a lower banded matrix 
           and its inverse if it exists.
*/


namespace details {
using namespace std;

template<typename Real>
class LowerBandedMatrix {
  
  using size_type = typename vector<Real>::size_type;

private:

  vector<size_type>     index;
  vector<vector<Real>>  band;
  size_type             num_bands;
  bool                  is_invertible;

public:   

  LowerBandedMatrix( vector<size_type> band_index, 
                     vector<vector<Real>> band_value ) : 
    index(band_index), band(band_value), 
    num_bands(index.size()), is_invertible(true) {
      
    auto it = find( index.begin(), index.end(), 0 );

    if( it == index.end() ) // no diagonal band
      is_invertible = false;
    else {
      for( auto e: band.at(0) ) {
        if( e == 0 ) {
          is_invertible = false;
          break;
        } // if( e == 0 )
      } // for( auto e: band[0] )
    } // else 
  } // end Constructor


  // Compute Ax += alpha*A*x for the range specified
  void apply( vector<Real>& Ax, 
              const vector<Real>& x,
              Real alpha,
              size_type begin, 
              size_type end ) {
 
    for( size_type row=begin; row<end; ++row ) {
      for( size_type i=0; i < num_bands; i++ ) {
        size_type col = row - index.at(i);
        if( row >= index.at(i) ) Ax.at(row) += alpha*x.at(col)*band.at(i).at(col-begin);
        else break;
      }
    }
  }

  void applyTranspose( vector<Real>& Ax,
                       const vector<Real>& x, 
                       Real alpha,
                       size_type begin,
                       size_type end ) {
    for( size_type row = begin; row<end; ++row ) {
      for( size_type i=0; i < num_bands; ++i ) {
        size_type col = row + index.at(i);
          if( row + index.at(i) < end ) Ax.at(row) += alpha*x.at(col)*band.at(i).at(row-begin);
          else break;
      }
    }
  }

  // Compute x += alpha*Ainv*(Ax) for the range specified
  void solve( vector<Real>& x, 
              const vector<Real>& Ax, 
              Real alpha,
              size_type begin, 
              size_type end ) {
    if( !is_invertible ) {
      throw logic_error("\nError: Cannot solve system. Matrix is singular.\n");
    }

    for( size_type row=begin; row<end; ++row ) {
      x.at(row) += alpha*Ax.at(row)/band.at(0).at(row);
      for( size_type i=1; i < num_bands; i++ ) {
        size_type col = row - index.at(i);
        if( row >= index.at(i) ) x.at(row) -= alpha*x.at(col)*band.at(i).at(col-begin);
        else break;
      }
    } 
  }


  void solveTranspose( vector<Real>& x, 
                       const vector<Real>& Ax, 
                       Real alpha,
                       size_type begin,
                       size_type end ) {
    if( !is_invertible ) {
      throw logic_error("\nError: Cannot solve system. Matrix is singular.\n");
    } 

    for( size_type row=begin; row<end; ++row ) {
      x[row] += alpha*Ax.at(row)/band.at(0).at(row);
      for( size_type i=1; i < num_bands; ++i ) {
        size_type col = row + index.at(i);
        if( row + index.at(i) < end ) x.at(row) -= alpha*x.at(col)*band.at(i).at(row-begin);                      
      }
    }
  }

  void print( ostream& os ) {
  
    int fieldWidth = 16;

    os << "Band indicies and values" << endl;
    for( auto i: index ) os << setw(16) << i;
    os << endl;
    
    size_type N = band.at(0).size();
    os << string(fieldWidth*3,'-') << endl;
    for( size_type l=0; l<N; ++l ) {
      for( size_type k=0; k<num_bands; ++k ) {
        if( l<band.at(k).size() )  os << setw(16) << band.at(k).at(l);
        else                    os << setw(16) << "*";
      }  
      os << endl;
    } 
  }

  const vector<Real>& at( size_type i ) const { return band.at(index.at(i)); }
  vector<Real>& at( size_type i ) { return band.at(index.at(i)); }

  const vector<Real>& operator[]( size_type i ) const { return band[index[i]]; }
  vector<Real>& operator[]( size_type i ) { return band[index[i]]; }

}; // class LowerBandedMatrix

} // namespace details

using details::LowerBandedMatrix;

#endif // LOWER_BANDED_MATRIX_HPP

