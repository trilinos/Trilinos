// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
            void join(Real &globalMin,
                      const Real &localMin) const {
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
            void join(int &globalFeasible,
                      const int &localFeasible) const {
                globalFeasible *= localFeasible;
            }

        }; // End struct Feasible

        //----------------------------------------------------------------------
        //
        // Project x onto the bounds
        template<class Real, class V, class CV>
        struct Project {
            typedef typename V::execution_space execution_space;
            V X_; // Optimization variable
            CV L_; // Lower bounds
            CV U_; // Upper bounds

            Project(V& X, const CV& L, const CV& U) : X_(X), L_(L), U_(U) {}

            KOKKOS_INLINE_FUNCTION
            void operator() (const int i) const {
                const int M = L_.extent(1);
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
        template<class Real, class V, class CV>
        struct PruneLowerActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            CV X_; // Optimization variable
            CV L_; // Lower bounds
            Real eps_;

            PruneLowerActive(V &Y, const CV &X, const CV &L,  Real eps) :
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
        template<class Real, class V, class CV>
        struct PruneUpperActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            CV X_; // Optimization variable
            CV U_; // Upper bounds
            Real eps_;

            PruneUpperActive(V &Y, const CV &X, const CV &U, Real eps) :
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
        template<class Real, class V, class CV>
        struct PruneActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            CV X_; // Optimization variable
            CV L_; // Lower bounds
            CV U_; // Upper bounds
            Real eps_;

            PruneActive(V &Y, const CV &X, const CV &L, const CV &U, Real eps) :
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
        template<class Real, class V, class CV>
        struct GradPruneLowerActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            CV G_; // Gradient
            CV X_; // Optimization variable
            CV L_; // Lower bounds

            Real xeps_, geps_;

            GradPruneLowerActive(V &Y, const CV &G, const CV &X, const CV &L,Real xeps, Real geps) :
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
        template<class Real, class V, class CV>
        struct GradPruneUpperActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            CV G_; // Gradient
            CV X_; // Optimization variable
            CV U_; // Upper bounds

            Real xeps_, geps_;

            GradPruneUpperActive(V &Y, const CV &G, const CV &X, const CV &U, Real xeps, Real geps) :
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
        template<class Real, class V, class CV>
        struct GradPruneActive {
            typedef typename V::execution_space execution_space;
            V Y_; // Variable to be pruned
            CV G_; // Gradient
            CV X_; // Optimization variable
            CV L_; // Lower bounds
            CV U_; // Upper bounds
            Real xeps_, geps_;

            GradPruneActive(V &Y, const CV &G, const CV &X, const CV &L, const CV &U, Real xeps, Real geps) :
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

    } // End namespace KokkosStructs


    //--------------------------------------------------------------------------


    template <class Real, class LO, class GO, class Node>
    class TpetraBoundConstraint : public BoundConstraint<Real> {


        typedef Tpetra::MultiVector<Real,LO,GO,Node> MV;
        typedef Tpetra::MultiVector<const Real,LO,GO,Node> CMV;
        typedef ROL::Ptr<MV> MVP;
        typedef TpetraMultiVector<Real,LO,GO,Node> TMV;
        typedef ROL::Ptr<TMV> TMVP;
        typedef typename MV::dual_view_type::t_dev ViewType;
        typedef typename CMV::dual_view_type::t_dev ConstViewType;

        private:
            int gblDim_;
            int lclDim_;
            ConstViewType l_;           // Kokkos view of Lower bounds
            ConstViewType u_;           // Kokkos view of Upper bounds
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
                l_(lp->getLocalViewDevice(Tpetra::Access::ReadOnly)),
                u_(up->getLocalViewDevice(Tpetra::Access::ReadOnly)),
                scale_(scale),
                comm_(lp->getMap()->getComm()) {

                KokkosStructs::MinGap<Real,ConstViewType> findmin(l_,u_);
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

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);

                KokkosStructs::Feasible<Real,ConstViewType> check(x_lcl, l_, u_);

                Kokkos::parallel_reduce(lclDim_,check,lclFeasible);

                Real gblFeasible;

                Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_MIN,lclFeasible,Teuchos::outArg(gblFeasible));

                return gblFeasible == 1 ? true : false;
            }

            void project( Vector<Real> &x ) {

                auto xp = getVector(x);

                ViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::Project<Real,ViewType,ConstViewType> proj(x_lcl,l_,u_);

                Kokkos::parallel_for(lclDim_,proj);
            }


            void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ViewType v_lcl = vp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::PruneLowerActive<Real,ViewType,ConstViewType> prune(v_lcl,x_lcl,l_,epsn);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ViewType v_lcl = vp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::PruneUpperActive<Real,ViewType,ConstViewType> prune(v_lcl,x_lcl,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) {
                auto xp = getVector(x);
                auto vp = getVector(v);

                Real epsn = std::min(scale_*eps,this->min_diff_);

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ViewType v_lcl = vp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::PruneActive<Real,ViewType,ConstViewType> prune(v_lcl,x_lcl,l_,u_,epsn);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
               Real epsn = std::min(scale_*xeps,this->min_diff_);

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ConstViewType g_lcl = gp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ViewType v_lcl = vp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::GradPruneLowerActive<Real,ViewType,ConstViewType> prune(v_lcl,g_lcl,x_lcl,l_,epsn,geps);

                Kokkos::parallel_for(lclDim_,prune);
            }

             void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*xeps,this->min_diff_);

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ConstViewType g_lcl = gp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ViewType v_lcl = vp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::GradPruneUpperActive<Real,ViewType,ConstViewType> prune(v_lcl,g_lcl,x_lcl,u_,epsn,geps);

                Kokkos::parallel_for(lclDim_,prune);
            }

            void pruneActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
                auto xp = getVector(x);
                auto gp = getVector(g);
                auto vp = getVector(v);
                Real epsn = std::min(scale_*xeps,this->min_diff_);

                ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ConstViewType g_lcl = gp->getLocalViewDevice(Tpetra::Access::ReadOnly);
                ViewType v_lcl = vp->getLocalViewDevice(Tpetra::Access::ReadWrite);

                KokkosStructs::GradPruneActive<Real,ViewType,ConstViewType> prune(v_lcl,g_lcl,x_lcl,l_,u_,epsn,geps);

                Kokkos::parallel_for(lclDim_,prune);
            }
      }; // End class TpetraBoundConstraint

} // End ROL Namespace

#endif
