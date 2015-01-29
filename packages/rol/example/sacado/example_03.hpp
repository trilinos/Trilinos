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

#include "ROL_StdVector.hpp"
#include "FiniteElement.hpp"
#include "LinearAlgebra.hpp"

using namespace ROL;


// Define some terms which make up the constraint 
class DiffU {
public:
    template<class ScalarT>
    ScalarT operator() (const ScalarT &u, const ScalarT &u_x, const ScalarT &z, const ScalarT &x) {
        return u_x;
    }    
};

class ExpU {
public:
    template<class ScalarT>
    ScalarT operator() (const ScalarT &u, const ScalarT &u_x, const ScalarT &z, const ScalarT &x) {
        return exp(u);
    }    
}; 

class MinusZ {
public:
    template<class ScalarT>
    ScalarT operator() (const ScalarT &u, const ScalarT &u_x, const ScalarT &z, const ScalarT &x) {
        return -z;
    }    
};



//! \brief Compute the constraint vector -u'' + exp(u) - z with natural boundary conditions
template<class Real>
class BoundaryValueProblem { 
    private:
        int n_; 
        VectorFunction<Real,DiffU>   diff_;       // -u'' term
        VectorFunction<Real,ExpU>    exp_;        // exp(u) term
        VectorFunction<Real,MinusZ>  minusz_;     // -z term    
    public:
        BoundaryValueProblem( Real xl, Real xr, Teuchos::RCP<NodalBasis<Real> > basisp  ) :
            n_(basisp->ni_), 
            diff_(xl,xr,true,basisp),
            exp_(xl,xr,false,basisp),
            minusz_(xl,xr,false,basisp) {
        }
        
        template<class ScalarT> 
        void value(Vector<ScalarT> &c, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {
             
            Teuchos::RCP<std::vector<ScalarT> > temp_rcp = Teuchos::rcp(new std::vector<ScalarT>(n_,0));
            StdVector<ScalarT> temp(temp_rcp);

            c.zero(); // This seems to be important. Getting weird results without this zeroing

            diff_.evaluate(u,z,c);       
            exp_.evaluate(u,z,temp);
            c.plus(temp); 
            minusz_.evaluate(u,z,temp);  
            c.plus(temp);                // c = -u" + exp(u) - z
        }
};




//! \brief Inherit and add method for applying the inverse partial constraint Jacobian and its adjoint 
template<class Real,template<class> class BoundaryValueProblem> 
class BVP_Constraint : public Sacado_EqualityConstraint_SimOpt<Real,BoundaryValueProblem> {
    private:
        Teuchos::LAPACK<int,Real> lapack_;
     
    public:

    BVP_Constraint(const BoundaryValueProblem<Real> &bvp) :
        Sacado_EqualityConstraint_SimOpt<Real,BoundaryValueProblem> (bvp){}         
        
    void applyInverseJacobian_1(Vector< Real > &ijv, const Vector< Real > &v,
                                const Vector< Real > &u, const Vector< Real > &z,  Real &tol) { 

        Teuchos::RCP<const std::vector<Real> > vp =
             (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();

        int n = vp->size();

        // Canonical vector
        Teuchos::RCP<std::vector<Real> > ep = Teuchos::rcp(new std::vector<Real>(n,0.0) ); 
        StdVector<Real> e(ep);

        // Jacobian applied to canonical vector
        Teuchos::RCP<std::vector<Real> > jep = Teuchos::rcp(new std::vector<Real>(n,0.0) ); 
        StdVector<Real> je(jep);

        Teuchos::RCP<std::vector<Real> > jacp = Teuchos::rcp(new std::vector<Real>(n*n,0.0));
        StdVector<Real> jac(jacp);
                
        // Form the column-stacked Jacobian 
        for(int i=0; i<n; ++i) {
            (*ep)[i] = 1.0;

            Sacado_EqualityConstraint_SimOpt<Real,BoundaryValueProblem>::applyJacobian_1(je,e,u,z,tol);
            std::copy(jep->begin(),jep->end(),jacp->begin()+n*i); 
            std::fill(ep->begin(),ep->end(),0.0);
        }               
   
        lusolve(lapack_,jac,v,ijv);
    }	

     void applyInverseAdjointJacobian_1(Vector< Real > &iajv, const Vector< Real > &v,
                                    const Vector< Real > &u, const Vector< Real > &z,  Real &tol) { 

        Teuchos::RCP<const std::vector<Real> > vp =
             (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();

        int n = vp->size();

        // Canonical vector
        Teuchos::RCP<std::vector<Real> > ep = Teuchos::rcp(new std::vector<Real>(n,0.0) ); 
        StdVector<Real> e(ep);

        // Adjoint Jacobian applied to canonical vector
        Teuchos::RCP<std::vector<Real> > ajep = Teuchos::rcp(new std::vector<Real>(n,0.0) ); 
        StdVector<Real> aje(ajep);

        Teuchos::RCP<std::vector<Real> > ajacp = Teuchos::rcp(new std::vector<Real>(n*n,0.0));
        StdVector<Real> ajac(ajacp);
                
        // Form the column-stacked Jacobian 
        for(int i=0; i<n; ++i) {
            (*ep)[i] = 1.0;

            Sacado_EqualityConstraint_SimOpt<Real,BoundaryValueProblem>::applyAdjointJacobian_1(aje,e,u,z,tol);
            std::copy(ajep->begin(),ajep->end(),ajacp->begin()+n*i); 
            std::fill(ep->begin(),ep->end(),0.0);
        }               
   
        lusolve(lapack_,ajac,v,iajv);
    }	


    void solve(Vector<Real> &u, const Vector<Real> &z, Real &tol) {
        Teuchos::RCP<std::vector<Real> > up = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());    

        int n = up->size();

        Teuchos::RCP<std::vector<Real> > cp = Teuchos::rcp(new std::vector<Real>(n,0.0) );
        StdVector<Real> c(cp);

        // \f$ -\delta u \f$
        Teuchos::RCP<std::vector<Real> > mdup = Teuchos::rcp(new std::vector<Real>(n,0.0) );
        StdVector<Real> mdu(mdup);

        Real res    = 1.0;
        Real restol = 1e-7;
        int MAXIT   = 20;
        int iter    = 0;

        // Newton's method
        while(iter<MAXIT && res>restol) {
            this->value(c,u,z,tol);
            this->applyInverseJacobian_1(mdu,c,u,z,tol); 
            mdu.scale(-1.0);
            u.plus(mdu);
            res = c.norm();
        }

        // Should throw exception if restol not met       
         
    }

};




//! \brief Objective to minimize \f[ \frac{1}{2} \|u-u_\text{targ}|\^2 + \frac{\gamma}{2}\|z\|^2 \f]
template<class Real>
class QuadraticTracking {
    private:
        Real gamma_;
        Teuchos::RCP<std::vector<Real> > u_targ_rcp_;
        Teuchos::RCP<NodalBasis<Real> > basisp_;
        int ni_;
        int nq_;
    public:
        QuadraticTracking( Real gamma, const Vector<Real> &u_targ, Teuchos::RCP<NodalBasis<Real> > basisp ) :
            gamma_(gamma), 
            u_targ_rcp_((Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(u_targ))).getVector() ), 
            basisp_(basisp), ni_(basisp_->ni_), nq_(basisp_->nq_) {}

        void update_gamma(Real gamma) { gamma_ = gamma; }
        void update_target(const Vector<Real> &u_targ) { 
            u_targ_rcp_ = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u_targ)).getVector()); } 

        //! \brief Compute \f[ \frac{1}{2} \|u-u_\text{targ}|\^2 + \frac{\gamma}{2}\|z\|^2 \f]
        template<class ScalarT>
        ScalarT value(const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {

            Teuchos::RCP<const std::vector<ScalarT> > up = 
               (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();
 
            Teuchos::RCP<const std::vector<ScalarT> > zp = 
               (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

            ScalarT err_norm = 0;
            ScalarT reg_norm = 0;

            // Numerically approximate the tracking and regularization norms
            for(int j=0; j<nq_; ++j) {
                ScalarT err_ptwise = 0;
                ScalarT reg_ptwise = 0;
                for(int i=0; i<ni_; ++i) {
                    err_ptwise += ((*up)[i]-(*u_targ_rcp_)[i])*(*basisp_->Lp_)[j+nq_*i];
                    reg_ptwise += (*zp)[i]*(*basisp_->Lp_)[j+nq_*i];
                }
                err_norm += 0.5*(*basisp_->wqp_)[j]*err_ptwise*err_ptwise;
                reg_norm += 0.5*(*basisp_->wqp_)[j]*reg_ptwise*reg_ptwise;
            }
            
            ScalarT J = err_norm + gamma_*reg_norm;

            return J;
        }        
};















