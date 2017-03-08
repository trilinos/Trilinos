// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
#ifndef MY_OPERATOR_HPP
#define MY_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"

//! Simple example of a user's defined Anasazi::Operator class.
/*! 
 * This is a simple, single processor example of user's defined
 * Anasazi::Operator-derived class. The class is templated with ScalarType;
 * possible choices are, for example, "float", "double", or
 * "complex<double>".
 *
 * This file can be easily extended to tackle more interesting cases.
 *
 * \author Oscar Chinallato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)
 *
 * \date Last modified on 01-Nov-05
 */
template <class ScalarType>
class MyOperator : public Anasazi::Operator<ScalarType>
{

public:
  
  /* Constructs a square matrix with \c NumRows rows and columns.
   * The matrix is tridiagonal, and the computational stencil is 
   * [-1, 2, -1]
   */
  MyOperator(const int NumRows) :
    NumRows_(NumRows)
  {
    l_ = -1.0;
    d_ =  2.0;
    u_ = -1.0;
  }

  // Constructor for tridiagonal matrix.
  MyOperator(const int NumRows, std::vector<ScalarType> ldu) :
    NumRows_(NumRows)
  {
    l_ = ldu[0];
    d_ = ldu[1];
    u_ = ldu[2];    
  }

  // Constructor for a diagonal matrix with variable entries.
  MyOperator(std::vector<ScalarType> diag) :
    NumRows_(diag.size())
  {
    diag_.resize(diag.size());
    for(unsigned int i=0; i<diag_.size(); ++i)
      diag_[i] = diag[i];
  }

  //! Dtor
  ~MyOperator()
  {}
  
  //! Applies the tridiagonal or diagonal matrix to a multivector.
  void Apply(const Anasazi::MultiVec<ScalarType>& X, 
                            Anasazi::MultiVec<ScalarType>& Y) const
  {
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    assert (MyX != 0);
    
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    assert (MyY != 0);
    
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetGlobalLength() == Y.GetGlobalLength());
   
    if (diag_.size() == 0)
    {
      // This is a tridiagonal matrix
      for (int v = 0 ; v < X.GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < X.GetGlobalLength() ; ++i)
        {
          if (i == 0)
          {
            (*MyY)[v][i] = (d_ * (*MyX)[v][i] + u_ * (*MyX)[v][i + 1]);
          }
          else if (i == X.GetGlobalLength() - 1)
          {
            (*MyY)[v][i] = (d_ * (*MyX)[v][i] + l_ * (*MyX)[v][i-1]);
          }
          else
          {
            (*MyY)[v][i] = (d_ * (*MyX)[v][i] + l_ * (*MyX)[v][i-1] + u_ * (*MyX)[v][i+1]);
          }
        }
      }
    } 
    else
    {
      // This is a diagonal matrix
      for (int v = 0 ; v < X.GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < X.GetGlobalLength() ; ++i)
        {
          (*MyY)[v][i] = diag_[i] * (*MyX)[v][i];
        }
      }      
    }
  }
  
private:
  //! Number of rows and columns
  int NumRows_;
  //! Elements on subdiagonal, diagonal, and superdiagonal.
  ScalarType l_, d_, u_;
  //! Elements on diagonal (for variable-diagonal case).
  std::vector<ScalarType> diag_;
};

#endif //MY_OPERATOR_HPP
