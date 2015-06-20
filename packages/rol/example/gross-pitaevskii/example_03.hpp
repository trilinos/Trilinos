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

/** \file   example_03.hpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence. In this version we implement the option to use 
            correct Sobolev inner products and Riesz mapping using nodal 
            Galerkin methods.         
                
    \details Minimize the one-dimensional Gross-Pitaevskii (GP) energy 
             functional
             \f[ J[\psi] = \int \frac{1}{2} |\nabla\psi|^2 + V(x)|\psi|^2 
                           +g|\psi|^4 \,\mathrm{d}x \f]
             Subject to the equality constraint that the particle density be
             normalized. 
             \f[ e(\psi) = \int |\psi|^2\,\mathrm{d}x - 1 = 0 \f]
             For simplicity, we will assume the wavefunction \f$\psi\f$ to 
             be real-valued, the potential function \f$ V(x)\geq 0\f$,
             the computational domain is the interval \f$[-1,1]\f$, and that
             \f$\psi(-1)=\psi(1)=0\f$. We descretize using the nodal Galerkin method
             \f[\psi(x)\approx \sum\limits_{k=1}^{n-1} \psi(x_k) \ell_k(x)\f].
               
    \author Greg von Winckel
    \date   Tue Dec 16 10:48:53 MST 2014
*/

#include <iostream>

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_CompositeStepSQP.hpp"
#include "ROL_Algorithm.hpp"

#include "numerics/NodalBasis.hpp"
#include "numerics/InnerProductMatrix.hpp"



using namespace ROL;

template <class Real, class Element=Real>
class OptStdVector;  // Optimization space.

template <class Real, class Element=Real>
class OptDualStdVector;  // Dual optimization space.

template <class Real, class Element=Real>
class ConStdVector;  // Constraint space.

template <class Real, class Element=Real>
class ConDualStdVector;  // Dual constraint space.

// Vector space definitions:


// Optimization space.
template <class Real, class Element>
class OptStdVector : public Vector<Real> {

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptDualStdVector<Real> >  dual_vec_;
bool useRiesz_; 
Teuchos::RCP<InnerProductMatrix<Real> > ipmat_;

public:

OptStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec, bool useRiesz, Teuchos::RCP<InnerProductMatrix<Real> > ipmat) : 
    std_vec_(std_vec), dual_vec_(Teuchos::null), useRiesz_(useRiesz), ipmat_(ipmat) {}

void plus( const Vector<Real> &x ) {
  OptStdVector &ex = Teuchos::dyn_cast<OptStdVector>(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  unsigned dimension = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}


//! Modify the dot product between primal variables to be \f$(u,v)=\int\limits_0^1 \dot u \dot v\,\mathrm{d}x \f$
Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  OptStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptStdVector<Real, Element> >(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  if(useRiesz_){ 
      val = this->ipmat_->inner(std_vec_,xvalptr);
  } 
  else{
      for (unsigned i=0; i<dimension; i++) {
          val += (*std_vec_)[i]*(*xvalptr)[i];
      }
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new OptStdVector( Teuchos::rcp( new std::vector<Element>(std_vec_->size()) ),useRiesz_,ipmat_ ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<OptStdVector> e = Teuchos::rcp( new OptStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)),useRiesz_, ipmat_ ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}


//! Modify the dual of vector u to be \f$\tilde u = -\ddot u\f$
const Vector<Real> & dual() const {

  dual_vec_ = Teuchos::rcp( new OptDualStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ), useRiesz_, ipmat_ ) );

  if(useRiesz_){
      this->ipmat_->apply(std_vec_, Teuchos::rcp_const_cast<std::vector<Real> >(dual_vec_->getVector()));
  }
  return *dual_vec_;
}

}; // class OptStdVector




// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public Vector<Real> {

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptStdVector<Real> >  dual_vec_;
bool useRiesz_;
Teuchos::RCP<InnerProductMatrix<Real> >ipmat_;

public:

OptDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec, bool useRiesz, Teuchos::RCP<InnerProductMatrix<Real> > ipmat) : 
    std_vec_(std_vec), dual_vec_(Teuchos::null), useRiesz_(useRiesz), ipmat_(ipmat) {}

void plus( const Vector<Real> &x ) {
  OptDualStdVector &ex = Teuchos::dyn_cast<OptDualStdVector>(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  unsigned dimension = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  OptDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<OptDualStdVector<Real, Element> >(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > kxvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  if(useRiesz_)
  {
      val = ipmat_->inv_inner(std_vec_,kxvalptr); 
  }
  else {
    for (unsigned i=0; i<dimension; i++) {
          val += (*std_vec_)[i]*(*kxvalptr)[i];
      }  

  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new OptDualStdVector( Teuchos::rcp( new std::vector<Element>(std_vec_->size()) ), useRiesz_, ipmat_ ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<OptDualStdVector> e = Teuchos::rcp( new OptDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)),useRiesz_, ipmat_ ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}

const Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new OptStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ), useRiesz_, ipmat_ ) );
  if(useRiesz_){
      this->ipmat_->solve(std_vec_, Teuchos::rcp_const_cast<std::vector<Real> >(dual_vec_->getVector()));
  }
  return *dual_vec_;
}

}; // class OptDualStdVector





// Constraint space.
template <class Real, class Element>
class ConStdVector : public Vector<Real> {

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<ConDualStdVector<Real> >  dual_vec_;

public:

ConStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const Vector<Real> &x ) {
  ConStdVector &ex = Teuchos::dyn_cast<ConStdVector>(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  unsigned dimension = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  ConStdVector<Real, Element> & ex = Teuchos::dyn_cast<ConStdVector<Real, Element> >(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();

   

  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new ConStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<ConStdVector> e = Teuchos::rcp( new ConStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)) ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}

const Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new ConDualStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class ConStdVector


// Dual constraint space.
template <class Real, class Element>
class ConDualStdVector : public Vector<Real> {
private:

Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<ConStdVector<Real> >  dual_vec_;

public:

ConDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const Vector<Real> &x ) {
  ConDualStdVector &ex = Teuchos::dyn_cast<ConDualStdVector>(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  unsigned dimension = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  ConDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<ConDualStdVector<Real, Element> >(const_cast <Vector<Real> &>(x));
  Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
  unsigned dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new ConDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<ConDualStdVector> e = Teuchos::rcp( new ConDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)) ) );
  (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
  return e;
}

int dimension() const {return std_vec_->size();}

const Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new ConStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class ConDualStdVector


/*** End of declaration of four vector spaces. ***/






/** Objective Function Class */
template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real> >
class Objective_GrossPitaevskii : public Objective<Real> {

    private:

        /** \var int ni_ Number of interior nodes  */ 

        const int  ni_;     
        const Real gnl_;
        Teuchos::RCP<NodalBasis<Real> > nb_;
         
        Teuchos::RCP<InnerProductMatrix<Real> > kinetic_;
        Teuchos::RCP<InnerProductMatrix<Real> > potential_;
        Teuchos::RCP<InnerProductMatrix<Real> > nonlinear_;

        void updateNonlinear(Teuchos::RCP<const std::vector<Real> > psip) {

             // This should be a member variable instead of reallocating each time
             Teuchos::RCP<std::vector<Real> > psi2q = 
                 Teuchos::rcp( new std::vector<Real> (nb_->nq_, 0.0) );

             nb_->lagrange_->dinterp(*psip,*psi2q);

             for(int i=0;i<nb_->nq_;++i) {
                 (*psi2q)[i] *= (*psi2q)[i];
             }
             nonlinear_->update(*psi2q);
        } 
        
    public: 

        Objective_GrossPitaevskii(const int ni, const Real gnl, 
                                  Teuchos::RCP<NodalBasis<Real> > nb,
                                  Teuchos::RCP<InnerProductMatrix<Real> > kinetic,
                                  Teuchos::RCP<InnerProductMatrix<Real> > potential,
                                  Teuchos::RCP<InnerProductMatrix<Real> > nonlinear):
                                  ni_(ni), gnl_(gnl), nb_(nb), kinetic_(kinetic),
                                  potential_(potential), nonlinear_(nonlinear) {}


   Real value( const Vector<Real> &psi, Real &tol ) {

        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        this->updateNonlinear(psip);

        // Compute objective function value
        Real J =  0.5*( kinetic_->inner(psip,psip) +  
                        potential_->inner(psip,psip) + 
                        gnl_*nonlinear_->inner(psip,psip) );
   
        return J;
    }


    void gradient( Vector<Real> &g, const Vector<Real> &psi, Real &tol ) {

        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        // Pointer to gradient vector 
        Teuchos::RCP<std::vector<Real> > gp = 
            Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<XDual>(g)).getVector());

        this->updateNonlinear(psip);

        kinetic_->apply(psip,gp);
        potential_->applyadd(psip,gp);     
        nonlinear_->applyaddtimes(psip,gp,2*gnl_);  
    }

    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol ) {
        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        // Pointer to direction vector 
        Teuchos::RCP<const std::vector<Real> > vp = (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(v))).getVector();

        // Pointer to action of Hessian on direction vector 
        Teuchos::RCP<std::vector<Real> > hvp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<XDual>(hv)).getVector());

        kinetic_->apply(vp,hvp);
        potential_->applyadd(vp,hvp);
        nonlinear_->applyaddtimes(vp,hvp,6*gnl_);

    }

};



/** Constraint class */
template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real>, class CPrim=StdVector<Real>, class CDual=StdVector<Real> >
class Normalization_Constraint : public EqualityConstraint<Real> {

    private:     
    Teuchos::RCP<InnerProductMatrix<Real> > mass_;

    public:
    Normalization_Constraint(Teuchos::RCP<InnerProductMatrix<Real> > mass) : 
        mass_(mass) {}          

    void value(Vector<Real> &c, const Vector<Real> &psi, Real &tol){

        // Pointer to constraint vector (only one element)
        Teuchos::RCP<std::vector<Real> > cp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<CPrim>(c)).getVector());

        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        (*cp)[0] = mass_->inner(psip,psip)-1.0;

    }

    void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

        // Pointer to action of Jacobian of constraint on direction vector (yields scalar)
        Teuchos::RCP<std::vector<Real> > jvp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<CPrim>(jv)).getVector());

        // Pointer to direction vector     
        Teuchos::RCP<const std::vector<Real> > vp = (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(v))).getVector();

        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        (*jvp)[0] = 2.0*mass_->inner(psip,vp);

    }

    void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

        // Pointer to action of adjoint of Jacobian of constraint on direction vector (yields vector)
        Teuchos::RCP<std::vector<Real> > ajvp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<XDual>(ajv)).getVector());

        // Pointer to direction vector     
        Teuchos::RCP<const std::vector<Real> > vp = (Teuchos::dyn_cast<CDual>(const_cast<Vector<Real> &>(v))).getVector();
 
        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        mass_->apply(psip,ajvp);
        for(unsigned i=0;i<psip->size();++i) {
            (*ajvp)[i] *= 2.0*(*vp)[0];
        }

    }

    void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, 
                             const Vector<Real> &psi, Real &tol){

        // The pointer to action of constraint Hessian in u,v inner product
        Teuchos::RCP<std::vector<Real> > ahuvp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<XDual>(ahuv)).getVector());

        // Pointer to direction vector u     
        Teuchos::RCP<const std::vector<Real> > up = (Teuchos::dyn_cast<CDual>(const_cast<Vector<Real> &>(u))).getVector();

        // Pointer to direction vector v     
        Teuchos::RCP<const std::vector<Real> > vp = (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(v))).getVector();

        // Pointer to opt vector 
        Teuchos::RCP<const std::vector<Real> > psip =
            (Teuchos::dyn_cast<XPrim>(const_cast<Vector<Real> &>(psi))).getVector();

        mass_->apply(vp,ahuvp);
        for(unsigned i=0;i<psip->size();++i) {
            (*ahuvp)[i] *= 2.0*(*up)[0];
        }
    }
     
};


