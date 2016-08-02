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

#ifndef PDE_PDEOBJECTIVE_HPP
#define PDE_PDEOBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_CompositeObjective_SimOpt.hpp"
#include "ROL_StdObjective.hpp"
#include "qoi.hpp"
#include "assembler.hpp"

template <class Real>
class PDE_Objective : public virtual ROL::Objective_SimOpt<Real> {
private:
  const std::vector<Teuchos::RCP<QoI<Real> > > qoi_vec_;
  const Teuchos::RCP<ROL::StdObjective<Real> > std_obj_;
  const Teuchos::RCP<Assembler<Real> > assembler_;

  std::vector<Teuchos::RCP<ROL::Objective_SimOpt<Real> > > obj_vec_;
  Teuchos::RCP<ROL::Objective_SimOpt<Real> > obj_;

  class IntegralObjective : public ROL::Objective_SimOpt<Real> {
    private:
      const Teuchos::RCP<QoI<Real> > qoi_;
      const Teuchos::RCP<Assembler<Real> > assembler_;
    public:
      IntegralObjective(const Teuchos::RCP<QoI<Real> > &qoi,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
        : qoi_(qoi), assembler_(assembler) {}

      Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
        return assembler_->assembleQoIValue(*up,*zp,*qoi_);
      }

      void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                      const ROL::Vector<Real> &z, Real &tol ) {
        Teuchos::RCP<Tpetra::MultiVector<> > gp =
          (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
        assembler_->assembleQoIGradient1(*up,*zp,*qoi_);
        gp->scale(static_cast<Real>(1),*(assembler_->getQoIGradient1()));
      }

      void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                      const ROL::Vector<Real> &z, Real &tol ) {
        Teuchos::RCP<Tpetra::MultiVector<> > gp =
          (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
        assembler_->assembleQoIGradient2(*up,*zp,*qoi_);
        gp->scale(static_cast<Real>(1),*(assembler_->getQoIGradient2()));
      }

      void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                 const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
        try {
          Teuchos::RCP<Tpetra::MultiVector<> > hvp =
            (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > vp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > up =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > zp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
          assembler_->assembleQoIHessVec11(*vp,*up,*zp,*qoi_);
          hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec11()));
        } catch (Exception::Zero &e) {
          hv.zero();
        }
      }
    
      void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                       const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
        try {
          Teuchos::RCP<Tpetra::MultiVector<> > hvp =
            (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > vp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > up =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > zp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
          assembler_->assembleQoIHessVec12(*vp,*up,*zp,*qoi_);
          hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec12()));
        } catch (Exception::Zero &e) {
          hv.zero();
        }
      }
    
      void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                       const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
        try {
          Teuchos::RCP<Tpetra::MultiVector<> > hvp =
            (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > vp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > up =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > zp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
          assembler_->assembleQoIHessVec21(*vp,*up,*zp,*qoi_);
          hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec21()));
        } catch (Exception::Zero &e) {
          hv.zero();
        }
      }
    
      void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                 const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
        try {
          Teuchos::RCP<Tpetra::MultiVector<> > hvp =
            (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > vp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > up =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
          Teuchos::RCP<const Tpetra::MultiVector<> > zp =
            (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
          assembler_->assembleQoIHessVec22(*vp,*up,*zp,*qoi_);
          hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec22()));
        } catch (Exception::Zero &e) {
          hv.zero();
        }
      }
  };

public:
  PDE_Objective(const std::vector<Teuchos::RCP<QoI<Real> > > &qoi_vec,
                const Teuchos::RCP<ROL::StdObjective<Real> > &std_obj,
                const Teuchos::RCP<Assembler<Real> > &assembler)
    : qoi_vec_(qoi_vec), std_obj_(std_obj), assembler_(assembler) {
    int size = qoi_vec_.size();
    obj_vec_.clear(); obj_vec_.resize(size,Teuchos::null);
    for (int i = 0; i < size; ++i) {
      obj_vec_[i] = Teuchos::rcp(new IntegralObjective(qoi_vec[i],assembler));
    }
    obj_ = Teuchos::rcp(new ROL::CompositeObjective_SimOpt<Real>(obj_vec_,std_obj_));
  }

  void update( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    obj_->update(u,z,flag,iter);
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    return obj_->value(u,z,tol);
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->gradient_1(g,u,z,tol);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->gradient_2(g,u,z,tol);
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->hessVec_11(hv,v,u,z,tol);
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->hessVec_12(hv,v,u,z,tol);
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->hessVec_21(hv,v,u,z,tol);
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->hessVec_22(hv,v,u,z,tol);
  }

}; // class PDE_Objective

#endif
