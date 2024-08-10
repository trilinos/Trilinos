// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
//
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TOP_ELE_WISE_BOUNDS_HELPERS_HPP
#define TOP_ELE_WISE_BOUNDS_HELPERS_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Thyra_VectorBase.hpp"
#include <cmath>

namespace RTOpPack {

  /** \brief Element-wise transformation operator for TOpEleWiseBound. */
  template<class Scalar>
    class TOpEleWiseBoundTransformation {
    public:
      void
      operator() (const Scalar &v0, const Scalar &v1, Scalar &z0) const {
        z0 = std::max (v0, std::min (v1, z0));
      }
    };

  /** \brief Element-wise transformation operator for TOpEleWisePrune. */
  template<class Scalar>
    class TOpEleWisePruneLowerTransformation {
    public:
    TOpEleWisePruneLowerTransformation(const Scalar& eps) : eps_(eps) {}
      void
      operator() (const Scalar &v0, const Scalar &v1, const Scalar &v2,
                  Scalar &z0) const {
        if ((v0 <= v1 + eps_) && (v2 > 0.0))
          z0 = 0;
      }
      void
      operator() (const Scalar &v0, const Scalar &v1, Scalar &z0) const {
        if (v0 <= v1 + eps_)
          z0 = 0;
      }

    private:
      const Scalar& eps_;
    };

  /** \brief Element-wise transformation operator for TOpEleWisePrune. */
  template<class Scalar>
    class TOpEleWisePruneUpperTransformation {
    public:
    TOpEleWisePruneUpperTransformation(const Scalar& eps) : eps_(eps) {}
      void
      operator() (const Scalar &v0, const Scalar &v1, const Scalar &v2,
                  Scalar &z0) const {
        if ((v0 >= v1 - eps_) && (v2 < 0.0))
          z0 = 0;
      }
      void
      operator() (const Scalar &v0, const Scalar &v1, Scalar &z0) const {
        if (v0 >= v1 - eps_)
          z0 = 0;
      }

    private:
      const Scalar& eps_;
    };


  /** \brief Element-wise product update transformation operator:
   * <tt>z0[i] *= min(v0[i],abs(z0[i]) * z0[i]/abs(z0[i]), i=0...n-1</tt>.
   */
  template<class Scalar>
    class TOpEleWiseBound : public TOp_2_1_Base<Scalar,
        TOpEleWiseBoundTransformation<Scalar> > {
    public:
      typedef TOp_2_1_Base<Scalar, TOpEleWiseBoundTransformation<Scalar> > base_t;
      /** \brief . */
      TOpEleWiseBound () :
          base_t (TOpEleWiseBoundTransformation<Scalar> ()) {
        this->setOpNameBase ("TOpElemWiseBound");
      }
    };



  /** \brief Element-wise product update transformation operator:
   * <tt>z0[i] *= min(v0[i],abs(z0[i]) * z0[i]/abs(z0[i]), i=0...n-1</tt>.
   */
  template<class Scalar>
    class TOpEleWisePruneLower_3_1 : public TOp_3_1_Base<Scalar,
        TOpEleWisePruneLowerTransformation<Scalar> > {
    public:
      typedef TOp_3_1_Base<Scalar, TOpEleWisePruneLowerTransformation<Scalar> > base_t;
      /** \brief . */
      TOpEleWisePruneLower_3_1 (const Scalar& eps) :
          base_t (TOpEleWisePruneLowerTransformation<Scalar> (eps)) {
        this->setOpNameBase ("TOpElemWisePruneLower_3_1");
      }
    };

  /** \brief Element-wise product update transformation operator:
   * <tt>z0[i] *= min(v0[i],abs(z0[i]) * z0[i]/abs(z0[i]), i=0...n-1</tt>.
   */
  template<class Scalar>
    class TOpEleWisePruneLower_2_1 : public TOp_2_1_Base<Scalar,
        TOpEleWisePruneLowerTransformation<Scalar> > {
    public:
      typedef TOp_2_1_Base<Scalar, TOpEleWisePruneLowerTransformation<Scalar> > base_t;
      /** \brief . */
      TOpEleWisePruneLower_2_1 (const Scalar& eps) :
          base_t (TOpEleWisePruneLowerTransformation<Scalar> (eps)) {
        this->setOpNameBase ("TOpElemWisePruneLower_2_1");
      }
    };

  /** \brief Element-wise product update transformation operator:
   * <tt>z0[i] *= min(v0[i],abs(z0[i]) * z0[i]/abs(z0[i]), i=0...n-1</tt>.
   */
  template<class Scalar>
    class TOpEleWisePruneUpper_3_1 : public TOp_3_1_Base<Scalar,
        TOpEleWisePruneUpperTransformation<Scalar> > {
    public:
      typedef TOp_3_1_Base<Scalar, TOpEleWisePruneUpperTransformation<Scalar> > base_t;
      /** \brief . */
      TOpEleWisePruneUpper_3_1 (const Scalar& eps) :
          base_t (TOpEleWisePruneUpperTransformation<Scalar> (eps)) {
        this->setOpNameBase ("TOpElemWisePruneUpper_3_1");
      }
    };

  /** \brief Element-wise product update transformation operator:
   * <tt>z0[i] *= min(v0[i],abs(z0[i]) * z0[i]/abs(z0[i]), i=0...n-1</tt>.
   */
  template<class Scalar>
    class TOpEleWisePruneUpper_2_1 : public TOp_2_1_Base<Scalar,
        TOpEleWisePruneUpperTransformation<Scalar> > {
    public:
      typedef TOp_2_1_Base<Scalar, TOpEleWisePruneUpperTransformation<Scalar> > base_t;
      /** \brief . */
      TOpEleWisePruneUpper_2_1 (const Scalar& eps) :
          base_t (TOpEleWisePruneUpperTransformation<Scalar> (eps)) {
        this->setOpNameBase ("TOpElemWisePruneUpper_2_1");
      }
    };



} // namespace RTOpPack

namespace Thyra {

  /** \brief Element-wise bound:
   * <tt>x(i) *= max(x_lo(i), min(x_up(i), x(i) ) ), i = 0...y->space()->dim()-1</tt>.
   *
   * \relates VectorBase
   */
  template<class Scalar>
    void
    ele_wise_bound (const ::Thyra::VectorBase<Scalar>& x_lo,
                       const ::Thyra::VectorBase<Scalar>& x_up,
                       const Teuchos::Ptr< ::Thyra::VectorBase<Scalar> > &x) {
      using Teuchos::tuple;
      using Teuchos::ptrInArg;
      using Teuchos::null;
      RTOpPack::TOpEleWiseBound<Scalar> ele_wise_bound_op;
      ::Thyra::applyOp<Scalar> (ele_wise_bound_op,
                                tuple (ptrInArg (x_lo), ptrInArg (x_up)), tuple (x),
                                Teuchos::null);
    }

  /** \brief Element-wise prune lower:
   * <tt> if ((x(i) <= x_lo(i) + eps_) && (g(i) > 0.0)), v(i) = 0, i = 0...y->space()->dim()-1</tt>.
   *
   * \relates VectorBase
   */
  template<class Scalar>
    void
    ele_wise_prune_lower (const ::Thyra::VectorBase<Scalar>& x,
                       const ::Thyra::VectorBase<Scalar>& x_lo,
                       const ::Thyra::VectorBase<Scalar>& g,
                       const Teuchos::Ptr< ::Thyra::VectorBase<Scalar> > &v,
                       const Scalar& eps) {
      using Teuchos::tuple;
      using Teuchos::ptrInArg;
      using Teuchos::null;
      RTOpPack::TOpEleWisePruneLower_3_1<Scalar> ele_wise_prune_op(eps);
      ::Thyra::applyOp<Scalar> (ele_wise_prune_op,
                                tuple (ptrInArg (x), ptrInArg (x_lo), ptrInArg (g)), tuple (v),
                                Teuchos::null);
    }

  /** \brief Element-wise prune lower:
   * <tt> if (x(i) <= x_lo(i) + eps_), v(i) = 0, i = 0...y->space()->dim()-1</tt>.
   *
   * \relates VectorBase
   */
  template<class Scalar>
    void
    ele_wise_prune_lower (const ::Thyra::VectorBase<Scalar>& x,
                       const ::Thyra::VectorBase<Scalar>& x_lo,
                       const Teuchos::Ptr< ::Thyra::VectorBase<Scalar> > &v,
                       const Scalar& eps) {
      using Teuchos::tuple;
      using Teuchos::ptrInArg;
      using Teuchos::null;
      RTOpPack::TOpEleWisePruneLower_2_1<Scalar> ele_wise_prune_op(eps);
      ::Thyra::applyOp<Scalar> (ele_wise_prune_op,
                                tuple (ptrInArg (x), ptrInArg (x_lo)), tuple (v),
                                Teuchos::null);
    }

  /** \brief Element-wise prune lower:
   * <tt> if ((x(i) >= x_up(i) - eps_) && (g(i) < 0.0)), v(i) = 0, i = 0...y->space()->dim()-1</tt>.
   *
   * \relates VectorBase
   */
  template<class Scalar>
    void
    ele_wise_prune_upper (const ::Thyra::VectorBase<Scalar>& x,
                       const ::Thyra::VectorBase<Scalar>& x_up,
                       const ::Thyra::VectorBase<Scalar>& g,
                       const Teuchos::Ptr< ::Thyra::VectorBase<Scalar> > &v,
                       const Scalar& eps) {
      using Teuchos::tuple;
      using Teuchos::ptrInArg;
      using Teuchos::null;
      RTOpPack::TOpEleWisePruneUpper_3_1<Scalar> ele_wise_prune_op(eps);
      ::Thyra::applyOp<Scalar> (ele_wise_prune_op,
                                tuple (ptrInArg (x), ptrInArg (x_up), ptrInArg (g)), tuple (v),
                                Teuchos::null);
    }

  /** \brief Element-wise prune lower:
   * <tt> if (x(i) >= x_up(i) - eps_), v(i) = 0, i = 0...y->space()->dim()-1</tt>.
   *
   * \relates VectorBase
   */

  template<class Scalar>
    void
    ele_wise_prune_upper (const ::Thyra::VectorBase<Scalar>& x,
                       const ::Thyra::VectorBase<Scalar>& x_up,
                       const Teuchos::Ptr< ::Thyra::VectorBase<Scalar> > &v,
                       const Scalar& eps) {
      using Teuchos::tuple;
      using Teuchos::ptrInArg;
      using Teuchos::null;
      RTOpPack::TOpEleWisePruneUpper_2_1<Scalar> ele_wise_prune_op(eps);
      ::Thyra::applyOp<Scalar> (ele_wise_prune_op,
                                tuple (ptrInArg (x), ptrInArg (x_up)), tuple (v),
                                Teuchos::null);
    }

}

#endif // RTOPPACK_TOP_ELE_WISE_BOUND_HPP
