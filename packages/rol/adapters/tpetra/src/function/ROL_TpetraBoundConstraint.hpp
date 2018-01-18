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
                const int M = L_.dimension_1();
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
                const int M = L_.dimension_1();
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
        // Project x onto the bounds  
        template<class Real, class V>
        struct Project {
            typedef typename V::execution_space execution_space;
            V X_; // Optimization variable
            V L_; // Lower bounds
            V U_; // Upper bounds

            Project(V& X, const V& L, const V& U) : X_(X), L_(L), U_(U) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const int M = L_.dimension_1();
                for(int j=0;j<M;++j) {
                    if( X_(i,j)<L_(i,j) ) {
                        X_(i,j) = L_(i,j); 
                    }
                    else if( X_(i,j)>U_(i,j) ) {
                        X_(i,j) = U_(i,j); 
                    }
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
                const int M = L_.dimension_1();
                for(int j=0;j<M;++j) {
                    if(X_(i,j)<=L_(i,j)+eps_) {
                        Y_(i,j) = 0.0;
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
                const int M = U_.dimension_1();
                for(int j=0;j<M;++j) {
                    if(X_(i,j)>=U_(i,j)-eps_) {
                        Y_(i,j) = 0.0;
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
                const int M = L_.dimension_1();
                for(int j=0;j<M;++j) {
                    if(X_(i,j)<=L_(i,j)+eps_) {
                        Y_(i,j) = 0.0;
                    } 
                    else if(X_(i,j)>=U_(i,j)-eps_) {
                        Y_(i,j) = 0.0;
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

            Real eps_;  

            GradPruneLowerActive(V &Y, const V &G, const V &X, const V &L,Real eps) : 
                Y_(Y), G_(G), X_(X), L_(L), eps_(eps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const int M = L_.dimension_1();
                for(int j=0;j<M;++j) {
                    if( (X_(i,j)<=L_(i,j)+eps_) && G_(i,j) > 0.0 ) {
                        Y_(i,j) = 0.0;
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

            Real eps_;  

            GradPruneUpperActive(V &Y, const V &G, const V &X, const V &U, Real eps) : 
                Y_(Y), G_(G), X_(X), U_(U), eps_(eps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const int M = U_.dimension_1();
                for(int j=0;j<M;++j) {
                    if( (X_(i,j)>=U_(i,j)-eps_) && G_(i,j) < 0.0 ) {
                        Y_(i,j) = 0.0;
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
            Real eps_;  

            GradPruneActive(V &Y, const V &G, const V &X, const V &L, const V &U, Real eps) : 
                Y_(Y), G_(G), X_(X), L_(L), U_(U), eps_(eps) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const int M = L_.dimension_1();
                for(int j=0;j<M;++j) {
                    if(( X_(i,j)<=L_(i,j)+eps_) && (G_(i,j)>0.0))  {
                        Y_(i,j) = 0.0;
                    } 
                    else if(( X_(i,j)>=U_(i,j)-eps_) && (G_(i,j)<0.0) ) {
                        Y_(i,j) = 0.0;
                    } 
                }
            }
        }; // End struct GradPruneActive

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
            Real min_diff_;
            Real scale_;
            ROL::Ptr<const Teuchos::Comm<int> > comm_; 

        ROL::Ptr<const MV> getVector( const ROL::Vector<Real>& x ) {
          return dynamic_cast<const TMV&>(x).getVector();
        }

        ROL::Ptr<MV> getVector( ROL::Vector<Real>& x ) {
          return dynamic_cast<TMV&>(x).getVector();
        }

        public:

            TpetraBoundConstraint(MVP lp, MVP up, Real scale = 1.0) :  
                gblDim_(lp->getGlobalLength()),
                lclDim_(lp->getLocalLength()),
                lp_(lp),
                up_(up),
                l_(lp->getDualView().d_view), 
                u_(up->getDualView().d_view), 
                scale_(scale),
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
                 
                ViewType x_lcl = xp->getDualView().d_view;
                
                KokkosStructs::Feasible<Real,ViewType> check(x_lcl, l_, u_);

                Kokkos::parallel_reduce(lclDim_,check,lclFeasible);

                Real gblFeasible;

                Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_MIN,lclFeasible,Teuchos::outArg(gblFeasible));

                return gblFeasible == 1 ? true : false;
            }            

            void project( Vector<Real> &x ) {

                auto xp = getVector(x);
                
                ViewType x_lcl = xp->getDualView().d_view;

                KokkosStructs::Project<Real,ViewType> proj(x_lcl,l_,u_);

                Kokkos::parallel_for(lclDim_,proj);
            }   
 

            void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getDualView().d_view;
                ViewType v_lcl = vp->getDualView().d_view;

                KokkosStructs::PruneLowerActive<Real,ViewType> prune(v_lcl,x_lcl,l_,epsn);

                Kokkos::parallel_for(lclDim_,prune);                
            }
              
            void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getDualView().d_view;
                ViewType v_lcl = vp->getDualView().d_view;

                KokkosStructs::PruneUpperActive<Real,ViewType> prune(v_lcl,x_lcl,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);                
            }
         
            void pruneActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getDualView().d_view;
                ViewType v_lcl = vp->getDualView().d_view;

                KokkosStructs::PruneActive<Real,ViewType> prune(v_lcl,x_lcl,l_,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);                
            }
         
            void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
               Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getDualView().d_view;
                ViewType g_lcl = gp->getDualView().d_view;
                ViewType v_lcl = vp->getDualView().d_view;

                KokkosStructs::GradPruneLowerActive<Real,ViewType> prune(v_lcl,g_lcl,x_lcl,l_,epsn);

                Kokkos::parallel_for(lclDim_,prune);                
            }
       
             void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getDualView().d_view;
                ViewType g_lcl = gp->getDualView().d_view;
                ViewType v_lcl = vp->getDualView().d_view;

                KokkosStructs::GradPruneUpperActive<Real,ViewType> prune(v_lcl,g_lcl,x_lcl,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);                
            }

            void pruneActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*eps,this->min_diff_);

                ViewType x_lcl = xp->getDualView().d_view;
                ViewType g_lcl = gp->getDualView().d_view;
                ViewType v_lcl = vp->getDualView().d_view;

                KokkosStructs::GradPruneActive<Real,ViewType> prune(v_lcl,g_lcl,x_lcl,l_,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);                
            }
      }; // End class TpetraBoundConstraint 

} // End ROL Namespace

#endif
