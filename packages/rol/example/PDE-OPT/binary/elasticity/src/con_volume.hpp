// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_QOI_VOLUME_HPP
#define PDEOPT_QOI_VOLUME_HPP

#include "ROL_StdConstraint.hpp"
#include "ROL_ParameterList.hpp"

template<class Real>
class Volume_TopoOpt : public ROL::StdConstraint<Real> {
private:
  Real w0_, vol_;
  std::vector<Real> w_;
  int N_, M_, T_;

public:
  Volume_TopoOpt(ROL::ParameterList &list) {
    w_  = ROL::getArrayFromStringParameter<Real>(list.sublist("Problem"), "Density");
    M_  = list.sublist("Problem").get("Number of Horizontal Cells",10);
    N_  = list.sublist("Problem").get("Number of Vertical Cells",20);
    T_  = w_.size();
    
    Real width  = list.sublist("Geometry").get("Width",2.0);
    Real height = list.sublist("Geometry").get("Height",1.0);
    vol_ = width*height/static_cast<Real>(M_*N_);

    Real mw(0);
    for (int i = 0; i < T_; ++i) {
      mw = std::max(mw,w_[i]);
    }
    w0_ = list.sublist("Problem").get("Maximum Weight Fraction", 0.5) * mw * width * height;
  }

  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    c[0] = -w0_;
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int k=0; k<T_; ++k) {
          c[0] += vol_ * w_[k] * x[i + M_*(j + N_*k)];
        }
      }
    }
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                const std::vector<Real> &x, Real &tol ) {
    jv[0] = static_cast<Real>(0);
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int k=0; k<T_; ++k) {
          jv[0] += vol_ * w_[k] * v[i + M_*(j + N_*k)];
        }
      }
    }
  }

   void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                        const std::vector<Real> &x, Real &tol ) {
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int k=0; k<T_; ++k) {
          ajv[i + M_*(j + N_*k)] = vol_ * w_[k] * v[0];
        }
      }
    }
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                      const std::vector<Real> &v, const std::vector<Real> &x,
                            Real &tol ) {
    int size = ahuv.size();
    ahuv.assign(size,static_cast<Real>(0));
  }
};

template<class Real>
class Selection_TopoOpt : public ROL::StdConstraint<Real> {
private:
  int N_, M_, T_;

public:
  Selection_TopoOpt(ROL::ParameterList &list) {
    std::vector<Real> w = ROL::getArrayFromStringParameter<Real>(list.sublist("Problem"), "Density");
    M_  = list.sublist("Problem").get("Number of Horizontal Cells",10);
    N_  = list.sublist("Problem").get("Number of Vertical Cells",20);
    T_  = w.size();
  }

  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        c[i + M_*j] = static_cast<Real>(-1);
        for (int k=0; k<T_; ++k) {
          c[i + M_*j] += x[i + M_*(j + N_*k)];
        }
      }
    }
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                const std::vector<Real> &x, Real &tol ) {
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        jv[i + M_*j] = static_cast<Real>(0);
        for (int k=0; k<T_; ++k) {
          jv[i + M_*j] += v[i + M_*(j + N_*k)];
        }
      }
    }
  }

   void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                        const std::vector<Real> &x, Real &tol ) {
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int k=0; k<T_; ++k) {
          ajv[i + M_*(j + N_*k)] = v[i + M_*j];
        }
      }
    }
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                      const std::vector<Real> &v, const std::vector<Real> &x,
                            Real &tol ) {
    int size = ahuv.size();
    ahuv.assign(size,static_cast<Real>(0));
  }
};

#endif
