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

/** \file
\brief  Contains definitions for Tpetra::MultiVector bound constraints.
\author Created by G. von Winckel
*/

#ifndef ROL_TPETRABOUNDCONSTRAINT_HPP
#define ROL_TPETRABOUNDCONSTRAINT_HPP

#include "Kokkos_Core.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_BoundConstraint.hpp"


namespace ROL {

    namespace KokkosStructs { // Parallel for and reduce functions

        //----------------------------------------------------------------------
        //
        // Find the minimum u_i-l_i
        template<class Real, class V>
        struct MinGap {
            typedef typename V::execution_space execution_space;
            V L_; // Lower bounds
            V U_; // Upper bounds

            MinGap(const V& L, const V& U) : L_(L), U_(U) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i, Real &min) const {
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    Real gap = U_(i,j)-L_(i,j);
                    if(gap<min) {
                        min = gap;
                    }
                }
            }

            KOKKOS_INLINE_FUNCTION
            void init(Real &min) const {
                min = U_(0,0)-L_(0,0);
            }

            KOKKOS_INLINE_FUNCTION
            void join(volatile Real &globalMin,
                      const volatile Real &localMin) const {
                if(localMin<globalMin) {
                    globalMin = localMin;
                }
            }
        }; // End struct MinGap

        //----------------------------------------------------------------------
        //
        // Determine if every l_i<=x_i<=u_i
        template<class Real, class V>
        struct Feasible {
            typedef typename V::execution_space execution_space;
            V X_; // Optimization variable
            V L_; // Lower bounds
            V U_; // Upper bounds

            Feasible(const V& X, const V& L, const V& U) : X_(X), L_(L), U_(U) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i, int &feasible) const {
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    if( (X_(i,j)<L_(i,j)) || (X_(i,j)>U_(i,j)) ) {
                        feasible = 0;
                    }
                }
            }

            KOKKOS_INLINE_FUNCTION
            void init(int &feasible) const {
                feasible = 1;
            }

            KOKKOS_INLINE_FUNCTION
            void join(volatile int &globalFeasible,
                      const volatile int &localFeasible) const {
                globalFeasible *= localFeasible;
            }

        }; // End struct Feasible

        //----------------------------------------------------------------------
        //
        // Helper functions
        template<class Real, class V>
        struct MinFunction {
            KOKKOS_INLINE_FUNCTION
            Real min(Real a,  Real b) const {
                return (a < b) ? a : b;
            }
        };   // End helper functions

        //----------------------------------------------------------------------
        //
        // Project x to the bounds
        template<class Real, class V>
        struct Project : MinFunction<Real, V> {
            typedef typename V::execution_space execution_space;
            V X_; // Optimization variable
            V L_; // Lower bounds
            V U_; // Upper bounds
            Real eps_, tol_, min_diff_;

            Project(V& X, const V& L, const V& U, Real eps = 0, Real tol = 0,
                                                  Real min_diff = 0) :
                X_(X), L_(L), U_(U), eps_(eps), tol_(tol), min_diff_(min_diff) {}

            KOKKOS_INLINE_FUNCTION
            Real max(Real a, Real b) const {
                return -this->min(-a, -b);
            }

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0), one(1);
                Real val1, val2;
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    val1 = L_(i,j);
                    if (eps_ != zero) {
                        val1 = ((L_(i,j) < -tol_) ? (one-eps_)*L_(i,j)
                             : ((L_(i,j) > +tol_) ? (one+eps_)*L_(i,j)
                             :   L_(i,j)+eps_));
                        val2 = L_(i,j)+eps_*min_diff_;
                        val1 = this->min(val1,val2);
                    }
                    X_(i,j) = max(X_(i,j), val1);

                    val1 = U_(i,j);
                    if (eps_ != zero) {
                        val1 = ((U_(i,j) < -tol_) ? (one+eps_)*U_(i,j)
                             : ((U_(i,j) > +tol_) ? (one-eps_)*U_(i,j)
                             :   U_(i,j)-eps_));
                        val2 = U_(i,j)-eps_*min_diff_;
                        val1 = max(val1,val2);
                    }
                    X_(i,j) = this->min(X_(i,j), val1);
                }
            }
        };   // End struct Project

        //----------------------------------------------------------------------
        //
        // Set variables to zero if they correspond to the lower active set
        template<class Real, class V>
        struct PruneLowerActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            V X_; // Optimization variable
            V L_; // Lower bounds
            Real eps_;

            PruneLowerActive(V &Y, const V &X, const V &L,  Real eps) :
                Y_(Y), X_(X), L_(L), eps_(eps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0);
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    if(X_(i,j)<=L_(i,j)+eps_) {
                        Y_(i,j) = zero;
                    }
                }
            }
        }; // End struct PruneLowerActive

        //----------------------------------------------------------------------
        //
        // Set variables to zero if they correspond to the upper active set
        template<class Real, class V>
        struct PruneUpperActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            V X_; // Optimization variable
            V U_; // Upper bounds
            Real eps_;

            PruneUpperActive(V &Y, const V &X, const V &U, Real eps) :
                Y_(Y), X_(X), U_(U), eps_(eps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0);
                const int M = U_.extent(1);
                for(int j=0;j<M;++j) {
                    if(X_(i,j)>=U_(i,j)-eps_) {
                        Y_(i,j) = zero;
                    }
                }
            }
        }; // End struct PruneUpperActive

        //----------------------------------------------------------------------
        //
        // Set variables to zero if they correspond to the active set
        template<class Real, class V>
        struct PruneActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            V X_; // Optimization variable
            V L_; // Lower bounds
            V U_; // Upper bounds
            Real eps_;

            PruneActive(V &Y, const V &X, const V &L, const V &U, Real eps) :
                Y_(Y), X_(X), L_(L), U_(U), eps_(eps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0);
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    if(X_(i,j)<=L_(i,j)+eps_) {
                        Y_(i,j) = zero;
                    }
                    else if(X_(i,j)>=U_(i,j)-eps_) {
                        Y_(i,j) = zero;
                    }
                }
            }
        }; // End struct PruneActive

        //----------------------------------------------------------------------
        //
        // Set variables to zero if they correspond to the lower active set and grad is positive
        template<class Real, class V>
        struct GradPruneLowerActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            V G_; // Gradient
            V X_; // Optimization variable
            V L_; // Lower bounds

            Real xeps_, geps_;

            GradPruneLowerActive(V &Y, const V &G, const V &X, const V &L,Real xeps, Real geps) :
                Y_(Y), G_(G), X_(X), L_(L), xeps_(xeps), geps_(geps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0);
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    if( (X_(i,j)<=L_(i,j)+xeps_) && (G_(i,j)>geps_) ) {
                        Y_(i,j) = zero;
                    }
                }
            }
        }; // End struct GradPruneLowerActive

        //----------------------------------------------------------------------
        //
        // Set variables to zero if they correspond to the upper active set and grad is negative
        template<class Real, class V>
        struct GradPruneUpperActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            V G_; // Gradient
            V X_; // Optimization variable
            V U_; // Upper bounds

            Real xeps_, geps_;

            GradPruneUpperActive(V &Y, const V &G, const V &X, const V &U, Real xeps, Real geps) :
                Y_(Y), G_(G), X_(X), U_(U), xeps_(xeps), geps_(geps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0);
                const int M = U_.extent(1);
                for(int j=0;j<M;++j) {
                    if( (X_(i,j)>=U_(i,j)-xeps_) && (G_(i,j)<-geps_) ) {
                        Y_(i,j) = zero;
                    }
                }
            }
        }; // End struct GradPruneUpperActive

        //----------------------------------------------------------------------
        //
        // Set variables to zero if they correspond to the active set
        template<class Real, class V>
        struct GradPruneActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            V G_; // Gradient
            V X_; // Optimization variable
            V L_; // Lower bounds
            V U_; // Upper bounds
            Real xeps_, geps_;

            GradPruneActive(V &Y, const V &G, const V &X, const V &L, const V &U, Real xeps, Real geps) :
                Y_(Y), G_(G), X_(X), L_(L), U_(U), xeps_(xeps), geps_(geps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const Real zero(0);
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    if(( X_(i,j)<=L_(i,j)+xeps_) && (G_(i,j)>geps_))  {
                        Y_(i,j) = zero;
                    }
                    else if(( X_(i,j)>=U_(i,j)-xeps_) && (G_(i,j)<-geps_) ) {
                        Y_(i,j) = zero;
                    }
                }
            }
        }; // End struct GradPruneActive

        //----------------------------------------------------------------------
        //
        // Build the scaling function for Coleman-Li.
        template<class Real, class V>
        struct BuildScalingFunction : MinFunction<Real,V> {
            typedef typename V::execution_space execution_space;
            V D_; // Scaling function
            V X_; // Optimization variable
            V G_; // Dual vector at which the scaling function is evaluated
            V L_; // Lower bounds
            V U_; // Upper bounds
            Real kappa_, zeta_;

            BuildScalingFunction(V &D, const V &X, const V &G,
                                       const V &L, const V &U,
                                       Real kappa, Real zeta) :
                D_(D), X_(X), G_(G), L_(L), U_(U), kappa_(kappa), zeta_(zeta)
                {}

            KOKKOS_INLINE_FUNCTION
            Real abs(Real a) const {
                return sgn(a)*a;
            }
            KOKKOS_INLINE_FUNCTION
            Real buildC(int i, int j) const {
                return this->min(zeta_*(U_(i,j) - L_(i,j)), kappa_);
            }
            KOKKOS_INLINE_FUNCTION
            Real sgn(Real a) const {
                return (a < 0) ? -1 : (a > 0 ? 1 : 0);
            }

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                Real c, lodiff, updiff;
                const int M = L_.extent(1);
                for(int j=0;j<M;++j) {
                    c = buildC(i,j);
                    lodiff = X_(i,j) - L_(i,j);
                    updiff = U_(i,j) - X_(i,j);
                    if (-G_(i,j) > lodiff) {
                        if (lodiff <= updiff) {
                            D_(i,j) = this->min(abs(G_(i,j)),c);
                            continue;
                        }
                    }
                    if (+G_(i,j) > updiff) {
                        if (updiff <= lodiff) {
                            D_(i,j) = this->min(abs(G_(i,j)),c);
                            continue;
                        }
                    }
                    D_(i,j) = this->min(lodiff, this->min(updiff, c));
                }
            }
        }; // End struct BuildScalingFunction

        //----------------------------------------------------------------------
        //
        // Apply the inverse of the scaling function for Coleman-Li.
        template<class Real, class V>
        struct ApplyInverseScalingFunction : BuildScalingFunction<Real,V> {
            V W_; // Primal vector being scaled

            ApplyInverseScalingFunction(V &D, const V &W, const V &X, const V &G,
                                              const V &L, const V &U,
                                              Real kappa, Real zeta) :
                BuildScalingFunction<Real,V>(D, X, G, L, U, kappa, zeta), W_(W) {};

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                BuildScalingFunction<Real,V>::operator()(i);
                const int M = this->L_.extent(1);
                for(int j=0;j<M;++j) {
                    this->D_(i,j) = this->W_(i,j)/this->D_(i,j);
                }
            }
        }; // End struct ApplyInverseScalingFunction

        //----------------------------------------------------------------------
        //
        // Apply the Jacobian of the scaling function for Coleman-Li.
        template<class Real, class V>
        struct ApplyScalingFunctionJacobian : BuildScalingFunction<Real,V> {
            V W_; // Primal vector being scaled

            ApplyScalingFunctionJacobian(V &D, const V &W, const V &X, const V &G,
                                               const V &L, const V &U,
                                               Real kappa, Real zeta) :
                BuildScalingFunction<Real,V>(D, X, G, L, U, kappa, zeta), W_(W) {};

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                BuildScalingFunction<Real,V>::operator()(i);
                const Real zero(0), one(1);
                Real indicator, d1prime;
                const int M = this->L_.extent(1);
                for(int j=0;j<M;++j) {
                    indicator = (this->D_(i,j) < this->buildC(i,j)) ? one : zero;
                    if (indicator == zero) {
                        this->D_(i,j) = zero;
                        continue;
                    }
                    d1prime = this->sgn(this->G_(i,j));
                    if (d1prime == zero) {
                        d1prime = one;
                        if (this->U_(i,j) - this->X_(i,j) < this->X_(i,j) - this->L_(i,j))
                            d1prime = -one;
                    }
                    this->D_(i,j) = d1prime*this->G_(i,j)*this->W_(i,j);
                }
            }
        }; // End struct ApplyScalingFunctionJacobian

    } // End namespace KokkosStructs


    //--------------------------------------------------------------------------


    template <class Real, class LO, class GO, class Node>
    class TpetraBoundConstraint : public BoundConstraint<Real> {

        typedef Tpetra::MultiVector<Real,LO,GO,Node> MV;
        typedef ROL::Ptr<MV> MVP;
        typedef ROL::Ptr<const MV> CMVP;
        typedef TpetraMultiVector<Real,LO,GO,Node> TMV;
        typedef ROL::Ptr<TMV> TMVP;
        typedef typename MV::dual_view_type::t_dev ViewType;

        private:
            int gblDim_;
            int lclDim_;
            MVP lp_;
            MVP up_;

            ViewType l_;           // Kokkos view of Lower bounds
            ViewType u_;           // Kokkos view of Upper bounds
            Real scale_;
            Real feasTol_;
            ROL::Ptr<const Teuchos::Comm<int> > comm_;
            Real min_diff_;

        ROL::Ptr<const MV> getVector( const ROL::Vector<Real>& x ) const {
          return dynamic_cast<const TMV&>(x).getVector();
        }

        ROL::Ptr<MV> getVector( ROL::Vector<Real>& x ) const {
          return dynamic_cast<TMV&>(x).getVector();
        }

        public:

            TpetraBoundConstraint(MVP lp, MVP up, Real scale = 1.0, 
                                  Real feasTol = std::sqrt(ROL_EPSILON<Real>())) :
                gblDim_(lp->getGlobalLength()),
                lclDim_(lp->getLocalLength()),
                lp_(lp),
                up_(up),
                l_(lp->getLocalViewDevice()),
                u_(up->getLocalViewDevice()),
                scale_(scale),
                feasTol_(feasTol),
                comm_(lp->getMap()->getComm()) {

                KokkosStructs::MinGap<Real,ViewType> findmin(l_,u_);
                Real lclMinGap = 0;

                // Reduce for this MPI process
                Kokkos::parallel_reduce(lclDim_,findmin,lclMinGap);

                Real gblMinGap;

                // Reduce over MPI processes
                Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_MIN,lclMinGap,Teuchos::outArg(gblMinGap));

                min_diff_ = 0.5*gblMinGap;
            }

            bool isFeasible( const Vector<Real> &x ) {
                auto xp = getVector(x);

                int lclFeasible = 1;

                ViewType x_lcl = xp->getLocalViewDevice();

                KokkosStructs::Feasible<Real,ViewType> check(x_lcl, l_, u_);

                Kokkos::parallel_reduce(lclDim_,check,lclFeasible);

                Real gblFeasible;

                Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_MIN,lclFeasible,Teuchos::outArg(gblFeasible));

                return gblFeasible == 1 ? true : false;
            }

            void project( Vector<Real> &x ) {
                auto xp = getVector(x);

                ViewType x_lcl = xp->getLocalViewDevice();

                KokkosStructs::Project<Real,ViewType> proj(x_lcl,l_,u_);

                Kokkos::parallel_for(lclDim_,proj);
            }

            void projectInterior( Vector<Real> &x ) {
                auto xp = getVector(x);

                ViewType x_lcl = xp->getLocalViewDevice();

                Real tol = 100.0*ROL_EPSILON<Real>();

                KokkosStructs::Project<Real,ViewType> projInterior(x_lcl,l_,u_,feasTol_,tol,min_diff_);

                Kokkos::parallel_for(lclDim_,projInterior);
            }

            void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getLocalViewDevice();
                ViewType v_lcl = vp->getLocalViewDevice();

                KokkosStructs::PruneLowerActive<Real,ViewType> prune(v_lcl,x_lcl,l_,epsn);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getLocalViewDevice();
                ViewType v_lcl = vp->getLocalViewDevice();

                KokkosStructs::PruneUpperActive<Real,ViewType> prune(v_lcl,x_lcl,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getLocalViewDevice();
                ViewType v_lcl = vp->getLocalViewDevice();

                KokkosStructs::PruneActive<Real,ViewType> prune(v_lcl,x_lcl,l_,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*xeps,this->min_diff_);

                ViewType x_lcl = xp->getLocalViewDevice();
                ViewType g_lcl = gp->getLocalViewDevice();
                ViewType v_lcl = vp->getLocalViewDevice();

                KokkosStructs::GradPruneLowerActive<Real,ViewType> prune(v_lcl,g_lcl,x_lcl,l_,epsn,geps);

                Kokkos::parallel_for(lclDim_,prune);
            }

             void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*xeps,this->min_diff_);

                ViewType x_lcl = xp->getLocalViewDevice();
                ViewType g_lcl = gp->getLocalViewDevice();
                ViewType v_lcl = vp->getLocalViewDevice();

                KokkosStructs::GradPruneUpperActive<Real,ViewType> prune(v_lcl,g_lcl,x_lcl,u_,epsn,geps);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*xeps,this->min_diff_);

                ViewType x_lcl = xp->getLocalViewDevice();
                ViewType g_lcl = gp->getLocalViewDevice();
                ViewType v_lcl = vp->getLocalViewDevice();

                KokkosStructs::GradPruneActive<Real,ViewType> prune(v_lcl,g_lcl,x_lcl,l_,u_,epsn,geps);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
                auto dvp = getVector(dv);
                auto vp  = getVector(v);
                auto xp  = getVector(x);
                auto gp  = getVector(g);
                Real kappa(1), zeta(0.5);

                ViewType dv_lcl = dvp->getLocalViewDevice();
                ViewType  v_lcl =  vp->getLocalViewDevice();
                ViewType  x_lcl =  xp->getLocalViewDevice();
                ViewType  g_lcl =  gp->getLocalViewDevice();

                KokkosStructs::ApplyInverseScalingFunction<Real,ViewType>
                    inverseScaling(dv_lcl,v_lcl,x_lcl,g_lcl,l_,u_,kappa,zeta);

                Kokkos::parallel_for(lclDim_,inverseScaling);
            }

            void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
                auto dvp = getVector(dv);
                auto vp  = getVector(v);
                auto xp  = getVector(x);
                auto gp  = getVector(g);
                Real kappa(1), zeta(0.5);

                ViewType dv_lcl = dvp->getLocalViewDevice();
                ViewType  v_lcl =  vp->getLocalViewDevice();
                ViewType  x_lcl =  xp->getLocalViewDevice();
                ViewType  g_lcl =  gp->getLocalViewDevice();

                KokkosStructs::ApplyScalingFunctionJacobian<Real,ViewType>
                    scalingJacobian(dv_lcl,v_lcl,x_lcl,g_lcl,l_,u_,kappa,zeta);

                Kokkos::parallel_for(lclDim_,scalingJacobian);
            }
      }; // End class TpetraBoundConstraint

} // End ROL Namespace

#endif
