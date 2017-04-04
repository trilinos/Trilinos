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

/*! \file  fieldhelper.hpp
    \brief
*/

#ifndef FIELDHELPER_HPP
#define FIELDHELPER_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_RCP.hpp"

template<class Real>
class FieldHelper {
  private:
  const int numFields_, numDofs_;
  const std::vector<int> numFieldDofs_;
  const std::vector<std::vector<int> > fieldPattern_;

  public:
  FieldHelper(const int numFields, const int numDofs,
              const std::vector<int> &numFieldDofs,
              const std::vector<std::vector<int> > &fieldPattern)
    : numFields_(numFields), numDofs_(numDofs),
      numFieldDofs_(numFieldDofs), fieldPattern_(fieldPattern) {}

  void splitFieldCoeff(std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & U,
                       const Teuchos::RCP<const Intrepid::FieldContainer<Real> >   & u_coeff) const {
    U.resize(numFields_);
    int  c = u_coeff->dimension(0);
    for (int i=0; i<numFields_; ++i) {
      U[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,numFieldDofs_[i]));
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          //U[i](j,k) = u_coeff(j,offset[i]+k);
          (*U[i])(j,k) = (*u_coeff)(j,fieldPattern_[i][k]);
        }
      }
    }
  }

  void combineFieldCoeff(Teuchos::RCP<Intrepid::FieldContainer<Real> >   & res,
                         const std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & R) const {
    int c = R[0]->dimension(0);  // number of cells
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, numDofs_));        
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          (*res)(j,fieldPattern_[i][k]) = (*R[i])(j,k);
        }
      }
    }
  }

  void combineFieldCoeff(Teuchos::RCP<Intrepid::FieldContainer<Real> >   & jac,
                         const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > & J) const {
    int c = J[0][0]->dimension(0);  // number of cells
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, numDofs_, numDofs_));        
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<numFields_; ++j) {
        for (int k=0; k<c; ++k) {
          for (int l=0; l<numFieldDofs_[i]; ++l) {
            for (int m=0; m<numFieldDofs_[j]; ++m) {
              (*jac)(k,fieldPattern_[i][l],fieldPattern_[j][m]) = (*J[i][j])(k,l,m);
            }
          }
        }
      }
    }
  }

  int numFields(void) const {
    return numFields_;
  }

};

#endif
