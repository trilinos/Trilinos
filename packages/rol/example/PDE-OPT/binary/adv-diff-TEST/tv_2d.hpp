// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  tv_2d.hpp
    \brief Provides the objective function definition for 2D total variation
	   with piecewise constant variables implemented using
	   TOOLS/Intrepid_HGRAD_C0_FEM.hpp.  This assumes the mesh cell
           ordering in the default rectangle mesh manager.
*/

#ifndef PDEOPT_QOI_ADV_DIFF_TV_2D_HPP
#define PDEOPT_QOI_ADV_DIFF_TV_2D_HPP

#include "../../TOOLS/pdevector.hpp"

template <class Real>
class Objective_TV_2D_C0 : public ROL::Objective<Real> {
private:
  Real volx_, voly_, vol_, eps_;
  int nx_, ny_;

public:
  Objective_TV_2D_C0(ROL::ParameterList & parlist) {
    Real width  = parlist.sublist("Geometry").get("Width",  2.0);
    Real height = parlist.sublist("Geometry").get("Height", 1.0);
    nx_ = parlist.sublist("Geometry").get("NX", 16);
    ny_ = parlist.sublist("Geometry").get("NY",  8);
    volx_ = width /static_cast<Real>(nx_);
    voly_ = height/static_cast<Real>(ny_);
    vol_  = volx_*voly_;
    eps_  = parlist.sublist("Problem").get("TV Smoothing Parameter",1e-3);
  }

  Real value(const ROL::Vector<Real> &x, Real &eps) {
    Teuchos::ArrayView<const Real> xd = (getConstField(x)->getData(0))();

    const Real half(0.5), two(2);
    Real sum(0);
    std::vector<Real> tmpx(nx_), tmpy(ny_);
    std::vector<std::vector<Real>> Dx(nx_+1,tmpy), Dy(ny_+1,tmpx);
    for (int i = 0; i < nx_+1; ++i) {
      for (int j = 0; j < ny_; ++j) {
        if (i==0) {
          Dx[i][j] = std::pow(-xd[i+j*nx_]/volx_,two);
        }
        else if (i > 0 && i < nx_) {
          Dx[i][j] = std::pow((xd[i+j*nx_]-xd[(i-1)+j*nx_])/volx_,two);
        }
        else {
          Dx[i][j] = std::pow(xd[(i-1)+j*nx_]/volx_,two);
        }
      }
    }
    for (int i = 0; i < ny_+1; ++i) {
      for (int j = 0; j < nx_; ++j) {
        if (i==0) {
          Dy[i][j] = std::pow(-xd[j+i*nx_]/voly_,two);
        }
        else if (i > 0 && i < ny_) {
          Dy[i][j] = std::pow((xd[j+i*nx_]-xd[j+(i-1)*nx_])/voly_,two);
        }
        else {
          Dy[i][j] = std::pow(xd[j+(i-1)*nx_]/voly_,two);
        }
      }
    }
    for (int i = 0; i < nx_+1; ++i) {
      for (int j = 0; j < ny_; ++j) {
        if (i==0) {
          sum += std::sqrt(Dx[i][j]+eps_);
        }
        else if (i > 0 && i < nx_) {
          sum += std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
        }
        else {
          sum += std::sqrt(half*Dx[i][j]+eps_);
        }
      }
    }
    for (int i = 0; i < ny_+1; ++i) {
      for (int j = 0; j < nx_; ++j) {
        if (i==0) {
          sum += std::sqrt(Dy[i][j]+eps_);
        }
        else if (i > 0 && i < ny_) {
          sum += std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
        }
        else {
          sum += std::sqrt(half*Dy[i][j]+eps_);
        }
      }
    }
    return vol_*sum;
  }

  void gradient(ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &eps) {
    g.zero();
    Teuchos::ArrayView<Real>       gd = (getField(g)->getDataNonConst(0))();
    Teuchos::ArrayView<const Real> xd = (getConstField(x)->getData(0))();

    const Real half(0.5), two(2), volx2 = volx_*volx_, voly2 = voly_*voly_;
    std::vector<Real> tmpx(nx_), tmpy(ny_);
    std::vector<std::vector<Real>> Dx(nx_+1,tmpy), Dy(ny_+1,tmpx);
    for (int i = 0; i < nx_+1; ++i) {
      for (int j = 0; j < ny_; ++j) {
        if (i==0) {
          Dx[i][j] = std::pow(-xd[i+j*nx_]/volx_,two);
        }
        else if (i > 0 && i < nx_) {
          Dx[i][j] = std::pow((xd[i+j*nx_]-xd[(i-1)+j*nx_])/volx_,two);
        }
        else {
          Dx[i][j] = std::pow(xd[(i-1)+j*nx_]/volx_,two);
        }
      }
    }
    for (int i = 0; i < ny_+1; ++i) {
      for (int j = 0; j < nx_; ++j) {
        if (i==0) {
          Dy[i][j] = std::pow(-xd[j+i*nx_]/voly_,two);
        }
        else if (i > 0 && i < ny_) {
          Dy[i][j] = std::pow((xd[j+i*nx_]-xd[j+(i-1)*nx_])/voly_,two);
        }
        else {
          Dy[i][j] = std::pow(xd[j+(i-1)*nx_]/voly_,two);
        }
      }
    }
    for (int i = 0; i < nx_+1; ++i) {
      for (int j = 0; j < ny_; ++j) {
        if (i==0) {
          Real cx = (vol_/volx2)/std::sqrt(Dx[i][j]+eps_);
          gd[i+j*nx_] += cx*xd[i+j*nx_];
        }
        else if (i > 0 && i < nx_-1) {
          Real cx = half*vol_/(volx2*std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_));
          gd[(i-1)+j*nx_] -= cx*(xd[i+j*nx_]-xd[(i-1)+j*nx_]);
          gd[i+j*nx_]     += cx*(two*xd[i+j*nx_]-xd[(i+1)+j*nx_]-xd[(i-1)+j*nx_]);
          gd[(i+1)+j*nx_] += cx*(xd[(i+1)+j*nx_]-xd[i+j*nx_]);
        }
        else if (i == nx_-1) {
          Real cx = half*vol_/(volx2*std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_));
          gd[(i-1)+j*nx_] -= cx*(xd[i+j*nx_]-xd[(i-1)+j*nx_]);
          gd[i+j*nx_]     += cx*(two*xd[i+j*nx_]-xd[(i-1)+j*nx_]);
        }
        else if (i==nx_) {
          Real cx = vol_*half/(volx2*std::sqrt(half*Dx[i][j]+eps_));
          gd[(i-1)+j*nx_] += cx*xd[(i-1)+j*nx_];
        }
      }
    }
    for (int i = 0; i < ny_+1; ++i) {
      for (int j = 0; j < nx_; ++j) {
        if (i==0) {
          Real cy = (vol_/voly2)/std::sqrt(Dy[i][j]+eps_);
          gd[j+i*nx_] += cy*xd[j+i*nx_];
        }
        else if (i > 0 && i < ny_-1) {
          Real cy = half*vol_/(voly2*std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_));
          gd[j+(i-1)*nx_] -= cy*(xd[j+i*nx_]-xd[j+(i-1)*nx_]);
          gd[j+i*nx_]     += cy*(two*xd[j+i*nx_]-xd[j+(i+1)*nx_]-xd[j+(i-1)*nx_]);
          gd[j+(i+1)*nx_] += cy*(xd[j+(i+1)*nx_]-xd[j+i*nx_]);
        }
        else if (i == ny_-1) {
          Real cy = half*vol_/(voly2*std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_));
          gd[j+(i-1)*nx_] -= cy*(xd[j+i*nx_]-xd[j+(i-1)*nx_]);
          gd[j+i*nx_]     += cy*(two*xd[j+i*nx_]-xd[j+(i-1)*nx_]);
        }
        else if (i==ny_) {
          Real cy = vol_*half/(voly2*std::sqrt(half*Dy[i][j]+eps_));
          gd[j+(i-1)*nx_] += cy*xd[j+(i-1)*nx_];
        }
      }
    }
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &eps) {
    hv.zero();
    Teuchos::ArrayView<Real>       hd = (getField(hv)->getDataNonConst(0))();
    Teuchos::ArrayView<const Real> vd = (getConstField(v)->getData(0))();
    Teuchos::ArrayView<const Real> xd = (getConstField(x)->getData(0))();

    const Real half(0.5), two(2), three(3), volx2 = volx_*volx_, voly2 = voly_*voly_;
    std::vector<Real> tmpx(nx_), tmpy(ny_);
    std::vector<std::vector<Real>> Dx(nx_+1,tmpy), Dy(ny_+1,tmpx);
    for (int i = 0; i < nx_+1; ++i) {
      for (int j = 0; j < ny_; ++j) {
        if (i==0) {
          Dx[i][j] = std::pow(-xd[i+j*nx_]/volx_,two);
        }
        else if (i > 0 && i < nx_) {
          Dx[i][j] = std::pow((xd[i+j*nx_]-xd[(i-1)+j*nx_])/volx_,two);
        }
        else {
          Dx[i][j] = std::pow(xd[(i-1)+j*nx_]/volx_,two);
        }
      }
    }
    for (int i = 0; i < ny_+1; ++i) {
      for (int j = 0; j < nx_; ++j) {
        if (i==0) {
          Dy[i][j] = std::pow(-xd[j+i*nx_]/voly_,two);
        }
        else if (i > 0 && i < ny_) {
          Dy[i][j] = std::pow((xd[j+i*nx_]-xd[j+(i-1)*nx_])/voly_,two);
        }
        else {
          Dy[i][j] = std::pow(xd[j+(i-1)*nx_]/voly_,two);
        }
      }
    }
    for (int i = 0; i < nx_+1; ++i) {
      for (int j = 0; j < ny_; ++j) {
        if (i==0) {
          Real cx1 = (vol_/volx2)/std::sqrt(Dx[i][j]+eps_);
          Real cx2 = (vol_/volx2)*(two/volx2)*(-half/std::pow(std::sqrt(Dx[i][j]+eps_),three));
          hd[i+j*nx_] += (cx1+cx2*std::pow(xd[i+j*nx_],two))*vd[i+j*nx_];
        }
        else if (i > 0 && i < nx_-1) {
          Real cx1 = (half*vol_/volx2)/std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
          Real cx2 = (half*vol_/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_),three));
          hd[(i-1)+j*nx_] += cx2*std::pow(xd[i+j*nx_]-xd[(i-1)+j*nx_],two)*vd[(i-1)+j*nx_]
                            -cx2*(xd[i+j*nx_]-xd[(i-1)+j*nx_])*(two*xd[i+j*nx_]-xd[(i+1)+j*nx_]-xd[(i-1)+j*nx_])*vd[i+j*nx_]
                            -cx2*(xd[i+j*nx_]-xd[(i-1)+j*nx_])*(xd[(i+1)+j*nx_]-xd[i+j*nx_])*vd[(i+1)+j*nx_]
                            -cx1*(vd[i+j*nx_]-vd[(i-1)+j*nx_]);
          hd[i+j*nx_]     += cx2*std::pow(two*xd[i+j*nx_]-xd[(i+1)+j*nx_]-xd[(i-1)+j*nx_],two)*vd[i+j*nx_]
                            -cx2*(two*xd[i+j*nx_]-xd[(i+1)+j*nx_]-xd[(i-1)+j*nx_])*(xd[i+j*nx_]-xd[(i-1)+j*nx_])*vd[(i-1)+j*nx_]
                            +cx2*(two*xd[i+j*nx_]-xd[(i+1)+j*nx_]-xd[(i-1)+j*nx_])*(xd[(i+1)+j*nx_]-xd[i+j*nx_])*vd[(i+1)+j*nx_]
                            +cx1*(two*vd[i+j*nx_]-vd[(i+1)+j*nx_]-vd[(i-1)+j*nx_]);
          hd[(i+1)+j*nx_] += cx2*std::pow(xd[(i+1)+j*nx_]-xd[i+j*nx_],two)*vd[(i+1)+j*nx_]
                            -cx2*(xd[(i+1)+j*nx_]-xd[i+j*nx_])*(xd[i+j*nx_]-xd[(i-1)+j*nx_])*vd[(i-1)+j*nx_]
                            +cx2*(xd[(i+1)+j*nx_]-xd[i+j*nx_])*(two*xd[i+j*nx_]-xd[(i+1)+j*nx_]-xd[(i-1)+j*nx_])*vd[i+j*nx_]
                            +cx1*(vd[(i+1)+j*nx_]-vd[i+j*nx_]);
        }
        else if (i == nx_-1) {
          Real cx1 = (half*vol_/volx2)/std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
          Real cx2 = (half*vol_/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_),three));
          hd[(i-1)+j*nx_] += cx2*std::pow(xd[i+j*nx_]-xd[(i-1)+j*nx_],two)*vd[(i-1)+j*nx_]
                            -cx2*(xd[i+j*nx_]-xd[(i-1)+j*nx_])*(two*xd[i+j*nx_]-xd[(i-1)+j*nx_])*vd[i+j*nx_]
                            -cx1*(vd[i+j*nx_]-vd[(i-1)+j*nx_]);
          hd[i+j*nx_]     += cx2*std::pow(two*xd[i+j*nx_]-xd[(i-1)+j*nx_],two)*vd[i+j*nx_]
                            -cx2*(two*xd[i+j*nx_]-xd[(i-1)+j*nx_])*(xd[i+j*nx_]-xd[(i-1)+j*nx_])*vd[(i-1)+j*nx_]
                            +cx1*(two*vd[i+j*nx_]-vd[(i-1)+j*nx_]);
        }
        else if (i==nx_) {
          Real cx1 = (vol_*half/volx2)/std::sqrt(half*Dx[i][j]+eps_);
          Real cx2 = (vol_*half/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*Dx[i][j]+eps_),three));
          hd[(i-1)+j*nx_] += (cx1+cx2*std::pow(xd[(i-1)+j*nx_],two))*vd[(i-1)+j*nx_];
        }
      }
    }
    for (int i = 0; i < ny_+1; ++i) {
      for (int j = 0; j < nx_; ++j) {
        if (i==0) {
          Real cy1 = (vol_/voly2)/std::sqrt(Dy[i][j]+eps_);
          Real cy2 = (vol_/voly2)*(two/voly2)*(-half/std::pow(std::sqrt(Dy[i][j]+eps_),three));
          hd[j+i*nx_] += (cy1+cy2*std::pow(xd[j+i*nx_],two))*vd[j+i*nx_];
        }
        else if (i > 0 && i < ny_-1) {
          Real cy1 = (half*vol_/voly2)/std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
          Real cy2 = (half*vol_/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_),three));
          hd[j+(i-1)*nx_] += cy2*std::pow(xd[j+i*nx_]-xd[j+(i-1)*nx_],two)*vd[j+(i-1)*nx_]
                            -cy2*(xd[j+i*nx_]-xd[j+(i-1)*nx_])*(two*xd[j+i*nx_]-xd[j+(i+1)*nx_]-xd[j+(i-1)*nx_])*vd[j+i*nx_]
                            -cy2*(xd[j+i*nx_]-xd[j+(i-1)*nx_])*(xd[j+(i+1)*nx_]-xd[j+i*nx_])*vd[j+(i+1)*nx_]
                            -cy1*(vd[j+i*nx_]-vd[j+(i-1)*nx_]);
          hd[j+i*nx_]     += cy2*std::pow(two*xd[j+i*nx_]-xd[j+(i+1)*nx_]-xd[j+(i-1)*nx_],two)*vd[j+i*nx_]
                            -cy2*(two*xd[j+i*nx_]-xd[j+(i+1)*nx_]-xd[j+(i-1)*nx_])*(xd[j+i*nx_]-xd[j+(i-1)*nx_])*vd[j+(i-1)*nx_]
                            +cy2*(two*xd[j+i*nx_]-xd[j+(i+1)*nx_]-xd[j+(i-1)*nx_])*(xd[j+(i+1)*nx_]-xd[j+i*nx_])*vd[j+(i+1)*nx_]
                            +cy1*(two*vd[j+i*nx_]-vd[j+(i+1)*nx_]-vd[j+(i-1)*nx_]);
          hd[j+(i+1)*nx_] += cy2*std::pow(xd[j+(i+1)*nx_]-xd[j+i*nx_],two)*vd[j+(i+1)*nx_]
                            -cy2*(xd[j+(i+1)*nx_]-xd[j+i*nx_])*(xd[j+i*nx_]-xd[j+(i-1)*nx_])*vd[j+(i-1)*nx_]
                            +cy2*(xd[j+(i+1)*nx_]-xd[j+i*nx_])*(two*xd[j+i*nx_]-xd[j+(i+1)*nx_]-xd[j+(i-1)*nx_])*vd[j+i*nx_]
                            +cy1*(vd[j+(i+1)*nx_]-vd[j+i*nx_]);
        }
        else if (i == ny_-1) {
          Real cy1 = (half*vol_/voly2)/std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
          Real cy2 = (half*vol_/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_),three));
          hd[j+(i-1)*nx_] += cy2*std::pow(xd[j+i*nx_]-xd[j+(i-1)*nx_],two)*vd[j+(i-1)*nx_]
                            -cy2*(xd[j+i*nx_]-xd[j+(i-1)*nx_])*(two*xd[j+i*nx_]-xd[j+(i-1)*nx_])*vd[j+i*nx_]
                            -cy1*(vd[j+i*nx_]-vd[j+(i-1)*nx_]);
          hd[j+i*nx_]     += cy2*std::pow(two*xd[j+i*nx_]-xd[j+(i-1)*nx_],two)*vd[j+i*nx_]
                            -cy2*(two*xd[j+i*nx_]-xd[j+(i-1)*nx_])*(xd[j+i*nx_]-xd[j+(i-1)*nx_])*vd[j+(i-1)*nx_]
                            +cy1*(two*vd[j+i*nx_]-vd[j+(i-1)*nx_]);
        }
        else if (i==ny_) {
          Real cy1 = (vol_*half/voly2)/std::sqrt(half*Dy[i][j]+eps_);
          Real cy2 = (vol_*half/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*Dy[i][j]+eps_),three));
          hd[j+(i-1)*nx_] += (cy1+cy2*std::pow(xd[j+(i-1)*nx_],two))*vd[j+(i-1)*nx_];
        }
      }
    }
  }

private:

  ROL::Ptr<const Tpetra::MultiVector<>> getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::TpetraMultiVector<Real>> xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
        if (xvec == ROL::nullPtr) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<Tpetra::MultiVector<>> getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::TpetraMultiVector<Real>> xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getField();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

}; // Objective_TV_2D_C0

#endif
