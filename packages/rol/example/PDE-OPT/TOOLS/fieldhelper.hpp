// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  fieldhelper.hpp
    \brief
*/

#ifndef FIELDHELPER_HPP
#define FIELDHELPER_HPP

#include "Intrepid_FieldContainer.hpp"
#include "ROL_Ptr.hpp"

namespace FieldUtils {

struct FieldInfo {
  const int numFields, numDofs;
  const std::vector<int> numFieldDofs;
  const std::vector<std::vector<int>> fieldPattern;

  FieldInfo(const int numFields_, const int numDofs_,
            const std::vector<int> &numFieldDofs_,
            const std::vector<std::vector<int>> &fieldPattern_)
    : numFields(numFields_), numDofs(numDofs_),
      numFieldDofs(numFieldDofs_), fieldPattern(fieldPattern_) {}
};

template<typename Real>
inline void splitFieldCoeff(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &U,
                            const ROL::Ptr<const Intrepid::FieldContainer<Real>>  &u_coeff,
                            const ROL::Ptr<const FieldInfo>                       &info) {
  const int numFields = info->numFields;
  U.resize(numFields);
  const int c = u_coeff->dimension(0);
  for (int i=0; i<numFields; ++i) {
    const int numFieldDofs = info->numFieldDofs[i];
    U[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,numFieldDofs);
    for (int j=0; j<c; ++j) {
      for (int k=0; k<numFieldDofs; ++k) {
        (*U[i])(j,k) = (*u_coeff)(j,info->fieldPattern[i][k]);
      }
    }
  }
}

template<typename Real>
inline void splitFieldCoeff(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &J,
                            const ROL::Ptr<const Intrepid::FieldContainer<Real>>               &jac,
                            const ROL::Ptr<const FieldInfo>                                    &rowInfo,
                            const ROL::Ptr<const FieldInfo>                                    &colInfo) {
  const int rowNumFields = rowInfo->numFields;
  const int colNumFields = colInfo->numFields;
  J.resize(rowNumFields);
  const int c = jac->dimension(0);
  for (int i=0; i<rowNumFields; ++i) {
    const int rowNumFieldDofs = rowInfo->numFieldDofs[i];
    J[i].resize(colNumFields,ROL::nullPtr);
    for (int j=0; j<colNumFields; ++j) {
      const int colNumFieldDofs = colInfo->numFieldDofs[j];
      J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,rowNumFieldDofs,colNumFieldDofs);
      for (int k=0; k<c; ++k) {
        for (int l=0; l<rowNumFieldDofs; ++l) {
          for (int m=0; m<colNumFieldDofs; ++m) {
            (*J[i][j])(k,l,m) = (*jac)(k,rowInfo->fieldPattern[i][l],colInfo->fieldPattern[j][m]);
          }
        }
      }
    }
  }
}

template<typename Real>
inline void combineFieldCoeff(ROL::Ptr<Intrepid::FieldContainer<Real>>                    &res,
                              const std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &R,
                              const ROL::Ptr<const FieldInfo>                             &info) {
  const int numFields = info->numFields;
  const int c = R[0]->dimension(0);  // number of cells
  res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, info->numDofs);
  for (int i=0; i<numFields; ++i) {
    const int numFieldDofs = info->numFieldDofs[i];
    for (int j=0; j<c; ++j) {
      for (int k=0; k<numFieldDofs; ++k) {
        (*res)(j,info->fieldPattern[i][k]) = (*R[i])(j,k);
      }
    }
  }
}

template<typename Real>
inline void combineFieldCoeff(ROL::Ptr<Intrepid::FieldContainer<Real>>                                 &jac,
                              const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &J,
                              const ROL::Ptr<const FieldInfo>                                          &rowInfo,
                              const ROL::Ptr<const FieldInfo>                                          &colInfo) {
  const int rowNumFields = rowInfo->numFields;
  const int colNumFields = colInfo->numFields;
  const int c = J[0][0]->dimension(0);  // number of cells
  jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, rowInfo->numDofs, colInfo->numDofs);
  for (int i=0; i<rowNumFields; ++i) {
    const int rowNumFieldDofs = rowInfo->numFieldDofs[i];
    for (int j=0; j<colNumFields; ++j) {
      const int colNumFieldDofs = colInfo->numFieldDofs[j];
      for (int k=0; k<c; ++k) {
        for (int l=0; l<rowNumFieldDofs; ++l) {
          for (int m=0; m<colNumFieldDofs; ++m) {
            (*jac)(k,rowInfo->fieldPattern[i][l],colInfo->fieldPattern[j][m]) = (*J[i][j])(k,l,m);
          }
        }
      }
    }
  }
}

}

template<class Real>
class FieldHelper {
  private:
  const int numFields_, numDofs_;
  const std::vector<int> numFieldDofs_;
  const std::vector<std::vector<int>> fieldPattern_;

  public:
  FieldHelper(const int numFields, const int numDofs,
              const std::vector<int> &numFieldDofs,
              const std::vector<std::vector<int>> &fieldPattern)
    : numFields_(numFields), numDofs_(numDofs),
      numFieldDofs_(numFieldDofs), fieldPattern_(fieldPattern) {}

  void splitFieldCoeff(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & U,
                       const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff) const {
    U.resize(numFields_);
    int  c = u_coeff->dimension(0);
    for (int i=0; i<numFields_; ++i) {
      U[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,numFieldDofs_[i]);
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          //U[i](j,k) = u_coeff(j,offset[i]+k);
          (*U[i])(j,k) = (*u_coeff)(j,fieldPattern_[i][k]);
        }
      }
    }
  }

  void splitFieldCoeff(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & J,
                       const ROL::Ptr<const Intrepid::FieldContainer<Real>> & jac) const {
    J.resize(numFields_);
    int  c = jac->dimension(0);
    for (int i=0; i<numFields_; ++i) {
      J[i].resize(numFields_,ROL::nullPtr);
      for (int j=0; j<numFields_; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,numFieldDofs_[i],numFieldDofs_[j]);
        for (int k=0; k<c; ++k) {
          for (int l=0; l<numFieldDofs_[i]; ++l) {
            for (int m=0; m<numFieldDofs_[j]; ++m) {
              (*J[i][j])(k,l,m) = (*jac)(k,fieldPattern_[i][l],fieldPattern_[j][m]);
            }
          }
        }
      }
    }
  }

  void combineFieldCoeff(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                         const std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & R) const {
    int c = R[0]->dimension(0);  // number of cells
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, numDofs_);        
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          (*res)(j,fieldPattern_[i][k]) = (*R[i])(j,k);
        }
      }
    }
  }

  void combineFieldCoeff(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                         const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & J) const {
    int c = J[0][0]->dimension(0);  // number of cells
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, numDofs_, numDofs_);        
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
