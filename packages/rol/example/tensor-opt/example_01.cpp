// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Contributed by Christoph Lohmann.
           Shows how to optimally limit tensor quantities; used in
           optimization-based flux correction schemes.
*/

////////////////////////////////////////////////////////////////////////////////////
// minimize F(x) = \sum_i (1 - x_i)^2 = (1 - x0)^2 + (1 - x1)^2 + (1 - x2)^2      //
// such that                                                                      //
//   0 <= x_i <= 1                                                                //
//   0 <=   A_i + F - lambda_i_min * I = B(  A_i,   F, - lambda_i_min)            //
//   0 <= - A_i - F + lambda_i_max * I = B(- A_i, - F,   lambda_i_max)            //
//   0 <=   A_j - F - lambda_j_min * I = B(  A_j, - F, - lambda_j_min)            //
//   0 <= - A_j + F + lambda_j_max * I = B(- A_j,   F,   lambda_j_max)            //
//                                                                                //
// where B(., ., .) is defined by                                                 //
//        B(A, F, lambda) := A + x .* F + lambda * I                              //
//       [a00, a01, a02;   [x0 * f0,       0,       0;   [lambda,      0,      0; //
//     =  a01, a11, a12; +        0, x1 * f1,       0; +       0, lambda,      0; //
//        a02, a12, a22]          0,       0, x2 * f2]         0,      0, lambda] //
////////////////////////////////////////////////////////////////////////////////////

#include <random>
#include <iostream>
#include <cassert>

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wunused-parameter"
//#pragma GCC diagnostic ignored "-Woverloaded-virtual"
//#pragma GCC diagnostic ignored "-Wsuggest-override"
//#pragma GCC diagnostic ignored "-Wterminate"
//#pragma GCC diagnostic ignored "-Wshadow"
//#pragma GCC diagnostic ignored "-Wundef"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "ROL_Vector.hpp"
#include "ROL_CArrayVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_StatusTest.hpp"

//#pragma GCC diagnostic pop

/// node i is the default one;
/// node j means we have to use -F due to F_{ij} = - F_{ji}
enum class Node
  {
   i,
   j,
  };

template <typename DT_, int dim_,
          typename std::enable_if<dim_ == 3, int>::type = 0>
class SemidefiniteProgramming
{
private:
  template <typename DT2_, typename DT1_>
  static auto my_cast(DT1_ & obj) -> decltype(dynamic_cast<DT2_>(obj))
  {
    return dynamic_cast<DT2_>(obj);
  }

private:
  /******************************************************************************/

  template <typename DT2_, int dim2_>
  class VectorWrapper : public ROL::Vector<DT2_>
  {
  private:
    bool _own;
    DT2_ * _pval;

  public:
    static constexpr int dim = dim2_;

  private:
    void _delete()
    {
      if (_own == true)
      {
        delete[] _pval;
      }
    }

  public:
    VectorWrapper() :
      _own(true),
      _pval(new DT2_[dim2_])
    {
    }

    VectorWrapper(const DT2_ & val_in) :
      VectorWrapper()
    {
      this->set(val_in);
    }

    VectorWrapper(const VectorWrapper & x) = delete;

    VectorWrapper(DT2_ * pval_in) :
      _own(false),
      _pval(pval_in)
    {
    }

    ~VectorWrapper()
    {
      _delete();
    }

    void wrap(DT2_ * pval_in)
    {
      _delete();
      _own = false;
      _pval = pval_in;
    }

    DT2_ * data()
    {
      return _pval;
    }

    DT2_ const * data() const
    {
      return _pval;
    }

    DT2_ & operator[](const int & i)
    {
      return _pval[i];
    }

    const DT2_ & operator[](const int & i) const
    {
      return _pval[i];
    }

    void plus(const ROL::Vector<DT2_> & x) override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dimension() != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vectors must have the same dimension." );
#endif

      auto & ex = my_cast<const VectorWrapper &>(x);

      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] += ex._pval[i];
      }
    }

    void scale(const DT2_ alpha) override
    {
      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] *= alpha;
      }
    }

    void scale(const DT2_ alpha, const VectorWrapper & x)
    {
      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] = alpha * x._pval[i];
      }
    }

    DT2_ dot(const ROL::Vector<DT2_> & x) const override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dimension() != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vectors must have the same dimension." );
#endif

      auto & ex = my_cast<const VectorWrapper &>(x);

      DT2_ tmp(0);
      for (int i(0); i < dim2_; ++i)
      {
        tmp += _pval[i] * ex._pval[i];
      }

      return tmp;
    }

    DT2_ norm() const override
    {
      DT2_ tmp(0);
      for (int i(0); i < dim2_; ++i)
      {
        tmp += _pval[i] * _pval[i];
      }

      return std::sqrt(tmp);
    }

    void axpy(const DT2_ alpha, const ROL::Vector<DT2_> & x) override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dimension() != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vectors must have the same dimension." );
#endif

      auto & ex = my_cast<const VectorWrapper &>(x);

      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] += alpha * ex._pval[i];
      }
    }

    ROL::Ptr<ROL::Vector<DT2_>> clone() const override
    {
      auto tmp = ROL::makePtr<VectorWrapper>();
      tmp->zero(); // TODO: necessary?
      return tmp;
    }

    void zero() override
    {
      this->set(DT2_(0));
    }

    void set(const DT2_ alpha)
    {
      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] = alpha;
      }
    }

    void set(const ROL::Vector<DT2_> & x) override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dimension() != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vectors must have the same dimension." );
#endif

      auto & ex = my_cast<const VectorWrapper &>(x);

      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] = ex._pval[i];
      }
    }

    int dimension() const override
    {
      return dim2_;
    }

    void applyUnary(const ROL::Elementwise::UnaryFunction<DT2_> & f) override
    {
      for (auto i(0); i < dim2_; ++i)
      {
        _pval[i] = f.apply(_pval[i]);
      }
    }

    void applyBinary(const ROL::Elementwise::BinaryFunction<DT2_> & f,
                     const ROL::Vector<DT2_> & x) override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dimension() != x.dimension(),
                                  std::invalid_argument,
                                 "Error: Vectors must have the same dimension." );
#endif

      const VectorWrapper & ex = my_cast<const VectorWrapper &>(x);
      for (auto i(0); i < dim2_; i++)
      {
        _pval[i] = f.apply(_pval[i], ex[i]);
      }
    }

    DT2_ reduce(const ROL::Elementwise::ReductionOp<DT2_> &r) const override
    {
      auto result = r.initialValue();
      for (auto i(0); i < dim2_; ++i)
      {
        r.reduce(_pval[i],result);
      }
      return result;
    }

    void print (std::ostream & outStream) const override
    {
      if (dim2_ == 0)
      {
        outStream << "[ ]" << std::endl;
      }
      else
      {
        outStream << "[" << _pval[0];
        for (int i(1); i < dim2_; ++i)
        {
          outStream << ", " << _pval[i];
        }
        outStream << "]" << std::endl;
      }
    }

    void random ()
    {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(1, 2);

      for (int i(0); i < dim2_; ++i)
      {
        _pval[i] = dis(gen);
      }
    }
  };

  /******************************************************************************/

  template <typename DT2_, int dim2_>
  class MyObjective : public ROL::Objective<DT2_>
  {
  private:
    const DT2_ _p;
  public:
    MyObjective(const DT2_ & p_in) :
      _p(p_in)
    {
      // at the moment method is restricted to _p = 2
      assert(p_in == DT2_(2));
    }

    DT2_ value (const ROL::Vector<DT2_> & x, DT2_ & /*tol*/) override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dim2_ != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vectors must have the same dimension." );
#endif

      auto & ex = my_cast<const VectorWrapper<DT2_, dim2_> &>(x);

      DT2_ res(0);

      for (int i(0); i < dim2_; ++i)
      {
        res += (DT2_(1) - ex[i]) * (DT2_(1) - ex[i]);
      }

      return res;
    }

    void gradient (ROL::Vector<DT2_> &       g,
                   const ROL::Vector<DT2_> & x,
                   DT2_ &                    /*tol*/) override
    {
#ifdef DEBUG
      ROL_TEST_FOR_EXCEPTION(dim2_ != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." );
      ROL_TEST_FOR_EXCEPTION(dim2_ != g.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." );
#endif

      auto & eg = my_cast<      VectorWrapper<DT2_, dim2_> &>(g);
      auto & ex = my_cast<const VectorWrapper<DT2_, dim2_> &>(x);

      for (int i(0); i < dim2_; ++i)
      {
        eg[i] = DT2_(2) * (ex[i] - DT2_(1));
      }
    }

    void hessVec (ROL::Vector<DT2_> &       hv,
                  const ROL::Vector<DT2_> & v,
                  const ROL::Vector<DT2_> & /*x*/,
                  DT2_ &                    /*tol*/) override
    {
#ifdef DEBUG
      /* ROL_TEST_FOR_EXCEPTION(dim2_ != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." ); */
      ROL_TEST_FOR_EXCEPTION(dim2_ != hv.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." );
      ROL_TEST_FOR_EXCEPTION(dim2_ != v.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." );
#endif

      auto & ehv = my_cast<      VectorWrapper<DT2_, dim2_> &>(hv);
      auto & ev  = my_cast<const VectorWrapper<DT2_, dim2_> &>(v);

      ehv.scale(DT2_(2), ev);
    }

    void invHessVec (ROL::Vector<DT2_> &       hv,
                     const ROL::Vector<DT2_> & v,
                     const ROL::Vector<DT2_> & /*x*/,
                     DT2_ &                    /*tol*/) override
    {
#ifdef DEBUG
      /* ROL_TEST_FOR_EXCEPTION(dim2_ != x.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." ); */
      ROL_TEST_FOR_EXCEPTION(dim2_ != hv.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." );
      ROL_TEST_FOR_EXCEPTION(dim2_ != v.dimension(),
                                 std::invalid_argument,
                                 "Error: Vector size is not correct." );
#endif

      auto & ehv = my_cast<      VectorWrapper<DT2_, dim2_> &>(hv);
      auto & ev  = my_cast<const VectorWrapper<DT2_, dim2_> &>(v);

      ehv.scale(DT2_(0.5), ev);
    }
  };

  /******************************************************************************/

  template <typename DT2_, int dim2_,
            typename std::enable_if<dim2_ == 3, int>::type = 0>
  class MyConstraint : public ROL::Constraint<DT2_>
  {
  private:
    VectorWrapper<DT2_, 6> _A;
    VectorWrapper<DT2_, 3> _F;
    static constexpr DT2_ _lambda = 0;

  public:
    void set_F(DT2_ * F_in)
    {
      _F.wrap(F_in);
    }

    void set_A(DT2_ * A_in)
    {
      _A.wrap(A_in);
    }

  public:
    static ROL::Ptr<VectorWrapper<DT2_, dim2_>> imul()
    {
      auto x = ROL::makePtr<VectorWrapper<DT2_, dim2_>>();
      x->zero();
      return x;
    }

    MyConstraint() :
      _A(nullptr),
      _F(nullptr)
    {
    }

    void value(ROL::Vector<DT2_> & c, const ROL::Vector<DT2_> & x, DT2_ & /*tol*/) override
    {
      auto & ex = my_cast<const VectorWrapper<DT2_, dim2_> &>(x);
      auto & ec = my_cast<VectorWrapper<DT2_, dim2_> &>(c);

      ec.zero();
      ec[0] = _A[0] + _A[3] + _A[5] + DT2_(3)*_lambda + _F[0]*ex[0] + _F[1]*ex[1] + _F[2]*ex[2];
      ec[1] = _A[0]*_A[3] + _A[0]*_A[5] + _A[3]*_A[5] + DT2_(2)*_A[0]*_lambda + DT2_(2)*_A[3]*_lambda
        + DT2_(2)*_A[5]*_lambda - _A[1]*_A[1] - _A[2]*_A[2] - _A[4]*_A[4]
        + DT2_(3)*_lambda*_lambda + _A[0]*_F[1]*ex[1] + _A[3]*_F[0]*ex[0] + _A[0]*_F[2]*ex[2]
        + _A[5]*_F[0]*ex[0] + _A[3]*_F[2]*ex[2] + _A[5]*_F[1]*ex[1] + DT2_(2)*_F[0]*_lambda*ex[0]
        + DT2_(2)*_F[1]*_lambda*ex[1] + DT2_(2)*_F[2]*_lambda*ex[2] + _F[0]*_F[1]*ex[0]*ex[1]
        + _F[0]*_F[2]*ex[0]*ex[2] + _F[1]*_F[2]*ex[1]*ex[2];
      ec[2] = _A[0]*_lambda*_lambda - _A[2]*_A[2]*_A[3] - _A[1]*_A[1]*_A[5]
        - _A[0]*_A[4]*_A[4] - _A[1]*_A[1]*_lambda - _A[2]*_A[2]*_lambda
        + _A[3]*_lambda*_lambda - _A[4]*_A[4]*_lambda + _A[5]*_lambda*_lambda
        + _lambda*_lambda*_lambda + _F[0]*_lambda*_lambda*ex[0] + _F[1]*_lambda*_lambda*ex[1]
        + _F[2]*_lambda*_lambda*ex[2] + DT2_(2)*_A[1]*_A[2]*_A[4] + _A[0]*_A[3]*_A[5]
        + _A[0]*_A[3]*_lambda + _A[0]*_A[5]*_lambda + _A[3]*_A[5]*_lambda
        - _A[2]*_A[2]*_F[1]*ex[1] - _A[4]*_A[4]*_F[0]*ex[0] - _A[1]*_A[1]*_F[2]*ex[2]
        + _A[0]*_F[1]*_lambda*ex[1] + _A[3]*_F[0]*_lambda*ex[0] + _A[0]*_F[2]*_lambda*ex[2]
        + _A[5]*_F[0]*_lambda*ex[0] + _A[3]*_F[2]*_lambda*ex[2] + _A[5]*_F[1]*_lambda*ex[1]
        + _A[0]*_A[3]*_F[2]*ex[2] + _A[0]*_A[5]*_F[1]*ex[1] + _A[3]*_A[5]*_F[0]*ex[0]
        + _A[0]*_F[1]*_F[2]*ex[1]*ex[2] + _A[3]*_F[0]*_F[2]*ex[0]*ex[2] + _A[5]*_F[0]*_F[1]*ex[0]*ex[1]
        + _F[0]*_F[1]*_lambda*ex[0]*ex[1] + _F[0]*_F[2]*_lambda*ex[0]*ex[2] + _F[1]*_F[2]*_lambda*ex[1]*ex[2]
        + _F[0]*_F[1]*_F[2]*ex[0]*ex[1]*ex[2];
       //ec.scale(1e2);
    }

    void applyJacobian(ROL::Vector<DT2_> &       jv,
                       const ROL::Vector<DT2_> & v,
                       const ROL::Vector<DT2_> & x,
                       DT2_ &                    /*tol*/) override
    {
      auto & ev  = my_cast<const VectorWrapper<DT2_, dim2_> &>(v);
      auto & ex  = my_cast<const VectorWrapper<DT2_, dim2_> &>(x);
      auto & ejv = my_cast<VectorWrapper<DT2_, dim2_> &>(jv);

      ejv.zero();
      ejv[0] += _F[0] * ev[0];
      ejv[0] += _F[1] * ev[1];
      ejv[0] += _F[2] * ev[2];

      ejv[1] += (_A[3]*_F[0] + _A[5]*_F[0] + DT2_(2)*_F[0]*_lambda + _F[0]*_F[1]*ex[1] + _F[0]*_F[2]*ex[2]) * ev[0];
      ejv[1] += (_A[0]*_F[1] + _A[5]*_F[1] + DT2_(2)*_F[1]*_lambda + _F[0]*_F[1]*ex[0] + _F[1]*_F[2]*ex[2]) * ev[1];
      ejv[1] += (_A[0]*_F[2] + _A[3]*_F[2] + DT2_(2)*_F[2]*_lambda + _F[0]*_F[2]*ex[0] + _F[1]*_F[2]*ex[1]) * ev[2];

      ejv[2] += (_F[0]*_lambda*_lambda - _A[4]*_A[4]*_F[0] + _A[3]*_A[5]*_F[0]
                 + _A[3]*_F[0]*_lambda + _A[5]*_F[0]*_lambda + _A[3]*_F[0]*_F[2]*ex[2] + _A[5]*_F[0]*_F[1]*ex[1]
                 + _F[0]*_F[1]*_lambda*ex[1] + _F[0]*_F[2]*_lambda*ex[2] + _F[0]*_F[1]*_F[2]*ex[1]*ex[2]) * ev[0];
      ejv[2] += (_F[1]*_lambda*_lambda - _A[2]*_A[2]*_F[1] + _A[0]*_A[5]*_F[1]
                 + _A[0]*_F[1]*_lambda + _A[5]*_F[1]*_lambda + _A[0]*_F[1]*_F[2]*ex[2] + _A[5]*_F[0]*_F[1]*ex[0]
                 + _F[0]*_F[1]*_lambda*ex[0] + _F[1]*_F[2]*_lambda*ex[2] + _F[0]*_F[1]*_F[2]*ex[0]*ex[2]) * ev[1];
      ejv[2] += (_F[2]*_lambda*_lambda - _A[1]*_A[1]*_F[2] + _A[0]*_A[3]*_F[2]
                 + _A[0]*_F[2]*_lambda + _A[3]*_F[2]*_lambda + _A[0]*_F[1]*_F[2]*ex[1] + _A[3]*_F[0]*_F[2]*ex[0]
                 + _F[0]*_F[2]*_lambda*ex[0] + _F[1]*_F[2]*_lambda*ex[1] + _F[0]*_F[1]*_F[2]*ex[0]*ex[1]) * ev[2];
      //ejv.scale(1e2);
    }

    void applyAdjointJacobian(ROL::Vector<DT2_> &       ajv,
                              const ROL::Vector<DT2_> & v,
                              const ROL::Vector<DT2_> & x,
                              DT2_ &                    /*tol*/) override
    {
      auto & ev   = my_cast<const VectorWrapper<DT2_, dim2_> &>(v);
      auto & ex   = my_cast<const VectorWrapper<DT2_, dim2_> &>(x);
      auto & eajv = my_cast<VectorWrapper<DT2_, dim2_> &>(ajv);

      eajv.zero();
      eajv[0] += _F[0] * ev[0];
      eajv[1] += _F[1] * ev[0];
      eajv[2] += _F[2] * ev[0];

      eajv[0] += (_A[3]*_F[0] + _A[5]*_F[0] + DT2_(2)*_F[0]*_lambda + _F[0]*_F[1]*ex[1] + _F[0]*_F[2]*ex[2]) * ev[1];
      eajv[1] += (_A[0]*_F[1] + _A[5]*_F[1] + DT2_(2)*_F[1]*_lambda + _F[0]*_F[1]*ex[0] + _F[1]*_F[2]*ex[2]) * ev[1];
      eajv[2] += (_A[0]*_F[2] + _A[3]*_F[2] + DT2_(2)*_F[2]*_lambda + _F[0]*_F[2]*ex[0] + _F[1]*_F[2]*ex[1]) * ev[1];

      eajv[0] += (_F[0]*_lambda*_lambda - _A[4]*_A[4]*_F[0] + _A[3]*_A[5]*_F[0]
                  + _A[3]*_F[0]*_lambda + _A[5]*_F[0]*_lambda + _A[3]*_F[0]*_F[2]*ex[2] + _A[5]*_F[0]*_F[1]*ex[1]
                  + _F[0]*_F[1]*_lambda*ex[1] + _F[0]*_F[2]*_lambda*ex[2] + _F[0]*_F[1]*_F[2]*ex[1]*ex[2]) * ev[2];
      eajv[1] += (_F[1]*_lambda*_lambda - _A[2]*_A[2]*_F[1] + _A[0]*_A[5]*_F[1]
                  + _A[0]*_F[1]*_lambda + _A[5]*_F[1]*_lambda + _A[0]*_F[1]*_F[2]*ex[2] + _A[5]*_F[0]*_F[1]*ex[0]
                  + _F[0]*_F[1]*_lambda*ex[0] + _F[1]*_F[2]*_lambda*ex[2] + _F[0]*_F[1]*_F[2]*ex[0]*ex[2]) * ev[2];
      eajv[2] += (_F[2]*_lambda*_lambda - _A[1]*_A[1]*_F[2] + _A[0]*_A[3]*_F[2]
                  + _A[0]*_F[2]*_lambda + _A[3]*_F[2]*_lambda + _A[0]*_F[1]*_F[2]*ex[1] + _A[3]*_F[0]*_F[2]*ex[0]
                  + _F[0]*_F[2]*_lambda*ex[0] + _F[1]*_F[2]*_lambda*ex[1] + _F[0]*_F[1]*_F[2]*ex[0]*ex[1]) * ev[2];
      //eajv.scale(1e2);
    }

    void applyAdjointJacobian(ROL::Vector<DT2_> &       ajv,
                              const ROL::Vector<DT2_> & v,
                              const ROL::Vector<DT2_> & x,
                              const ROL::Vector<DT2_> & /*dualv*/,
                              DT2_ &                    tol) override
    {
      applyAdjointJacobian(ajv, v, x, tol);
    }

    void applyAdjointHessian(ROL::Vector<DT2_> & ahuv,
                             const ROL::Vector<DT2_> & u,
                             const ROL::Vector<DT2_> & v,
                             const ROL::Vector<DT2_> & x,
                             DT2_ & /*tol*/) override
    {
      auto & ev    = my_cast<const VectorWrapper<DT2_, dim2_> &>(v);
      auto & eu    = my_cast<const VectorWrapper<DT2_, dim2_> &>(u);
      auto & ex    = my_cast<const VectorWrapper<DT2_, dim2_> &>(x);
      auto & eahuv = my_cast<VectorWrapper<DT2_, dim2_> &>(ahuv);

      eahuv.zero();
      eahuv[0] += (DT2_(0) * ev[0] + _F[0]*_F[1] * ev[1] + _F[0]*_F[2] * ev[2]) * eu[1];
      eahuv[1] += (_F[0]*_F[1] * ev[0] + DT2_(0) * ev[1] + _F[1]*_F[2] * ev[2]) * eu[1];
      eahuv[2] += (_F[0]*_F[2] * ev[0] + _F[1]*_F[2] * ev[1] + DT2_(0) * ev[2]) * eu[1];

      eahuv[0] += (DT2_(0) * ev[0]
                   + (_A[5]*_F[0]*_F[1] + _F[0]*_F[1]*_lambda + _F[0]*_F[1]*_F[2]*ex[2]) * ev[1]
                   + (_A[3]*_F[0]*_F[2] + _F[0]*_F[2]*_lambda + _F[0]*_F[1]*_F[2]*ex[1]) * ev[2]) * eu[2];
      eahuv[1] += ((_A[5]*_F[0]*_F[1] + _F[0]*_F[1]*_lambda + _F[0]*_F[1]*_F[2]*ex[2]) * ev[0]
                   + DT2_(0) * ev[1]
                   + (_A[0]*_F[1]*_F[2] + _F[1]*_F[2]*_lambda + _F[0]*_F[1]*_F[2]*ex[0]) * ev[2]) * eu[2];
      eahuv[2] += ((_A[3]*_F[0]*_F[2] + _F[0]*_F[2]*_lambda + _F[0]*_F[1]*_F[2]*ex[1]) * ev[0]
                   + (_A[0]*_F[1]*_F[2] + _F[1]*_F[2]*_lambda + _F[0]*_F[1]*_F[2]*ex[0]) * ev[1]
                   + DT2_(0) * ev[2]) * eu[2];
      //eahuv.scale(1e2);
    }
  };

  /******************************************************************************/

private:
  const std::string                                    _parfile;
  const ROL::Ptr<ROL::ParameterList>           _parlist;

  const ROL::Ptr<VectorWrapper<DT_, dim_>>         _lower,  _upper;
  const ROL::Ptr<VectorWrapper<DT_, dim_>>         _x;

  const ROL::Ptr<ROL::BoundConstraint<DT_>>        _bnd;

  const ROL::Ptr<ROL::BoundConstraint<DT_>>        _ibnd_local;

  std::vector<ROL::Ptr<ROL::Constraint<DT_>>>      _icon;
  std::vector<ROL::Ptr<ROL::Vector<DT_>>>          _imul;
  std::vector<ROL::Ptr<ROL::BoundConstraint<DT_>>> _ibnd;

  const ROL::Ptr<MyObjective<DT_, dim_>>           _obj;
  ROL::Ptr<ROL::OptimizationProblem<DT_>>          _problem;
  ROL::Ptr<ROL::OptimizationSolver<DT_>>           _solver;

  DT_ _A_i_up[6];
  DT_ _A_i_lo[6];
  DT_ _A_j_up[6];
  DT_ _A_j_lo[6];
  DT_ _F_p[3];
  DT_ _F_n[3];

public:
  /// Constructor
  SemidefiniteProgramming() :
    _parfile("example_01.xml"),
    _parlist(   Teuchos::rcp( new ROL::ParameterList() )),
    _lower(     ROL::makePtr<VectorWrapper<DT_, dim_>>(DT_(0))),
    _upper(     ROL::makePtr<VectorWrapper<DT_, dim_>>(DT_(1))),
    _x(         ROL::makePtr<VectorWrapper<DT_, dim_>>()),
    _bnd(       ROL::makePtr<ROL::Bounds<DT_>>(_lower, _upper)),
    _ibnd_local(ROL::makePtr<ROL::Bounds<DT_>>(* _lower, /*isLower =*/ true)),
    _obj(       ROL::makePtr<MyObjective<DT_, dim_>>(DT_(2)))
  {
    Teuchos::updateParametersFromXmlFile(_parfile, _parlist.ptr());

    for (int i(0); i < 4; ++i)
    {
      _icon.push_back(ROL::makePtr<MyConstraint<DT_, dim_>>());
      _imul.push_back(MyConstraint<DT_, dim_>::imul());
      _ibnd.push_back(_ibnd_local);
    }

    my_cast<MyConstraint<DT_, dim_> &>(* _icon[0]).set_A(_A_i_up);
    my_cast<MyConstraint<DT_, dim_> &>(* _icon[0]).set_F(_F_n);

    my_cast<MyConstraint<DT_, dim_> &>(* _icon[1]).set_A(_A_i_lo);
    my_cast<MyConstraint<DT_, dim_> &>(* _icon[1]).set_F(_F_p);

    my_cast<MyConstraint<DT_, dim_> &>(* _icon[2]).set_A(_A_j_up);
    my_cast<MyConstraint<DT_, dim_> &>(* _icon[2]).set_F(_F_p);

    my_cast<MyConstraint<DT_, dim_> &>(* _icon[3]).set_A(_A_j_lo);
    my_cast<MyConstraint<DT_, dim_> &>(* _icon[3]).set_F(_F_n);

    _problem = ROL::makePtr<ROL::OptimizationProblem<DT_>>(_obj, _x, _bnd, _icon, _imul, _ibnd);
    _solver = ROL::makePtr<ROL::OptimizationSolver<DT_>>(* _problem, * _parlist);
    _x->zero();
  }

  template <Node node_,
            typename std::enable_if<node_ == Node::i, int>::type = 0>
  void set_node(const DT_ A_in[6], const DT_ & lambda_min, const DT_ & lambda_max)
  {
    _A_i_up[0] = lambda_max - A_in[0];
    _A_i_up[1] =            - A_in[1];
    _A_i_up[2] =            - A_in[2];
    _A_i_up[3] = lambda_max - A_in[3];
    _A_i_up[4] =            - A_in[4];
    _A_i_up[5] = lambda_max - A_in[5];

    _A_i_lo[0] = A_in[0] - lambda_min;
    _A_i_lo[1] = A_in[1];
    _A_i_lo[2] = A_in[2];
    _A_i_lo[3] = A_in[3] - lambda_min;
    _A_i_lo[4] = A_in[4];
    _A_i_lo[5] = A_in[5] - lambda_min;
  }

  template <Node node_,
            typename std::enable_if<node_ == Node::j, int>::type = 0>
  void set_node(const DT_ A_in[6], const DT_ & lambda_min, const DT_ & lambda_max)
  {
    _A_j_up[0] = lambda_max - A_in[0];
    _A_j_up[1] =            - A_in[1];
    _A_j_up[2] =            - A_in[2];
    _A_j_up[3] = lambda_max - A_in[3];
    _A_j_up[4] =            - A_in[4];
    _A_j_up[5] = lambda_max - A_in[5];

    _A_j_lo[0] = A_in[0] - lambda_min;
    _A_j_lo[1] = A_in[1];
    _A_j_lo[2] = A_in[2];
    _A_j_lo[3] = A_in[3] - lambda_min;
    _A_j_lo[4] = A_in[4];
    _A_j_lo[5] = A_in[5] - lambda_min;
  }

  void set_flux(DT_ F_in[3])
  {
    for (int i(0); i < 3; ++i)
    {
      _F_p[i] = F_in[i];
    }
    for (int i(0); i < 3; ++i)
    {
      _F_n[i] = - F_in[i];
    }
  }

  DT_ * solve(DT_ x[3], std::ostream & outStream = std::cout)
  {
    _x->wrap(x);
    for (auto& it : _imul) it->zero();
    _solver->reset();
    _problem->reset();
    _solver->solve(outStream);

    return _x->data();
  }

  void check(std::ostream & outStream = std::cout)
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(DT_(- 0.1), DT_(0.1));

    DT_ A_i[] = {2.0, 0.0, 0.0
                 ,    1.0, 0.0
                 ,         0.0};
    DT_ A_j[] = {1.0, 0.0, 0.0
                 ,    2.0, 0.0
                 ,         0.0};

    DT_ F[] = {-2.0, -2.0, 0.5};
    DT_ lambda_i_up(3);
    DT_ lambda_i_lo(0);
    DT_ lambda_j_up(3);
    DT_ lambda_j_lo(0);

    set_node<Node::i>(A_i, lambda_i_lo, lambda_i_up);
    set_node<Node::j>(A_j, lambda_j_lo, lambda_j_up);
    set_flux(F);

    _problem->check(outStream);
  }

  void checkConstraints(DT_ sol[3])
  {
    int errorFlag(0);
    DT_ tol(1e-10);

    _x->wrap(sol);

    if (_bnd->isFeasible(* _x) == false)
    {
      ++errorFlag;
    }

    auto c = _x->clone();

    for (unsigned int i(0); i < _ibnd.size(); ++i)
    {
      DT_ eps(0);
      _icon[i]->value(* c, my_cast<ROL::Vector<DT_> &>(* _x), eps);

      for (int j(0); j < dim_; ++j)
      {
        if (my_cast<VectorWrapper<DT_, dim_> &>(* c)[j] < - tol)
        {
          ++errorFlag;
        }
      }
    }
    assert(errorFlag == 0);
  }
}; // end SemidefiniteProgramming



// ================================================================================
// ================================================================================

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  using DataType = double;
  constexpr int dim = 3;

  try
  {
    SemidefiniteProgramming<DataType, dim> sp;
    sp.check();

    {
      DataType A_i[] = {2.0, 0.0, 0.0
                        ,    1.0, 0.0
                        ,         0.0};
      DataType A_j[] = {1.0, 0.0, 0.0
                        ,    2.0, 0.0
                        ,         0.0};

      DataType F[] = {-2.0, -2.0, 0.5};
      DataType lambda_i_up(3);
      DataType lambda_i_lo(0);
      DataType lambda_j_up(3);
      DataType lambda_j_lo(0);

      sp.set_node<Node::i>(A_i, lambda_i_lo, lambda_i_up);
      sp.set_node<Node::j>(A_j, lambda_j_lo, lambda_j_up);
      sp.set_flux(F);
      // test exact solution
      DataType x[] = {1.0, 0.5, 0.0};
      sp.checkConstraints(x);
      sp.solve(x);
      std::cout << std::setprecision(16) << "x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
      // start from zero solution
      DataType y[] = {0.0, 0.0, 0.0};
      sp.solve(y);
      std::cout << "y = [" << y[0] << ", " << y[1] << ", " << y[2] << "]" << std::endl;
      // solve one more time
      DataType z[] = {0.0, 0.0, 0.0};
      sp.solve(z);
      std::cout << "z = [" << z[0] << ", " << z[1] << ", " << z[2] << "]" << std::endl;
      // perform checks
      ROL::CArrayVector<DataType> xx(&x[0],dim), yy(&y[0],dim), zz(&z[0],dim);
      xx.axpy(static_cast<DataType>(-1), yy);
      if (xx.norm() > std::sqrt(ROL::ROL_EPSILON<DataType>())) {
        *outStream << "\n\nxx.norm() = " << xx.norm() << "\n"; 
        errorFlag = 1000;
      }
      yy.axpy(static_cast<DataType>(-1), zz);
      if (yy.norm() > ROL::ROL_EPSILON<DataType>()) {
        *outStream << "\n\nyy.norm() = " << yy.norm() << "\n"; 
        errorFlag = 1000;
      }
    }
  }
  catch (std::logic_error& err)
  {
    std::cout << err.what() << std::endl;
    return 1;
  };

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
