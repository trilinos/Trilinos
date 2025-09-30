// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_TPETRA_1DFEM_FUNCTORS_HPP
#define NOX_TPETRA_1DFEM_FUNCTORS_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"

template <class ViewType, class GO>
struct IntegralOperatorReductionFunctor
{
  typedef typename ViewType::traits::non_const_value_type Scalar;

  const ViewType u_;
  const GO myMinGID_;
  const GO globalRow_;

  KOKKOS_INLINE_FUNCTION
  IntegralOperatorReductionFunctor(const ViewType& u,
                                   const GO& myMinGID,
                                   const std::size_t& globalRow) :
    u_(u),
    myMinGID_(myMinGID),
    globalRow_(globalRow)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const std::size_t i, Scalar& localSum) const
  {
    GO myGlobalRow = myMinGID_ + i;
    localSum += (globalRow_ + 0.5)*u_(i, 0)/static_cast<Scalar>(globalRow_ + myGlobalRow + 1);
  }

};

//==================================================================

template<class TpetraVectorType>
struct IntegralOperatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type Scalar;
  typedef typename TpetraVectorType::local_ordinal_type LO;
  typedef typename TpetraVectorType::global_ordinal_type GO;
  typedef typename TpetraVectorType::dual_view_type::t_dev ViewType;
  typedef typename TpetraVectorType::dual_view_type::t_dev::const_type ConstViewType;
  typedef typename TpetraVectorType::execution_space ExecSpace;
  typedef Kokkos::TeamPolicy<ExecSpace> TeamPolicy;
  typedef typename TeamPolicy::member_type MemberType;

  ConstViewType uView_;
  ViewType resultView_;
  const GO myMinGID_;
  const GO procMinGID_;
  const std::size_t myLength_;

  IntegralOperatorFunctor(const TpetraVectorType& u,
                          TpetraVectorType& result,
                          const GO& myMinGID,
                          const GO& procMinGID) :
    uView_(u.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly)),
    resultView_(result.template getLocalView<ExecSpace>(Tpetra::Access::ReadWrite)),
    myMinGID_(myMinGID),
    procMinGID_(procMinGID),
    myLength_(u.getLocalLength())
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const MemberType& member) const
  {
    int localRow = member.league_rank();
    const GO globalRow = procMinGID_ + localRow;

    Scalar localSum = 0;
    IntegralOperatorReductionFunctor<ConstViewType, GO> functor(uView_, myMinGID_, globalRow);
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(member, myLength_), functor, localSum);

    if (member.team_rank() == 0)
      resultView_(localRow, 0) = localSum;
  }

};

//==================================================================

template <class TpetraVectorType>
struct ResidualEvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type Scalar;
  typedef typename TpetraVectorType::execution_space ExecSpace;
  typedef typename TpetraVectorType::dual_view_type::t_dev ViewType;
  typedef typename TpetraVectorType::dual_view_type::t_dev::const_type ConstViewType;

  ViewType fView_;
  ConstViewType uView_;
  ConstViewType integralOpView_;

  ResidualEvaluatorFunctor(TpetraVectorType& f,
                           const TpetraVectorType& u,
                           const TpetraVectorType& integralOp) :
    fView_(f.template getLocalView<ExecSpace>(Tpetra::Access::ReadWrite)),
    uView_(u.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly)),
    integralOpView_(integralOp.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly))
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const std::size_t row) const
  {
    Scalar one = Kokkos::ArithTraits<Scalar>::one();
    fView_(row, 0) = one/(one - integralOpView_(row,0)) - uView_(row, 0);
  }

};

//==================================================================

template <class TpetraVectorType>
struct JacobianEvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type Scalar;
  typedef typename TpetraVectorType::dual_view_type::t_dev ViewType;
  typedef typename TpetraVectorType::dual_view_type::t_dev::const_type ConstViewType;
  typedef typename TpetraVectorType::execution_space ExecSpace;

  ViewType yView_;
  ConstViewType xView_;
  ConstViewType integralOpView_;
  ConstViewType integralOpXView_;
  const Scalar alpha_;
  const Scalar beta_;
  const Scalar omega_;
  const Tpetra::global_size_t globalLength_;

  JacobianEvaluatorFunctor(TpetraVectorType& y,
                           const TpetraVectorType& x,
                           const TpetraVectorType& integralOp,
                           const TpetraVectorType& integralOpX,
                           const Scalar& alpha,
                           const Scalar& beta,
                           const Scalar& omega) :
    yView_(y.template getLocalView<ExecSpace>(Tpetra::Access::ReadWrite)),
    xView_(x.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly)),
    integralOpView_(integralOp.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly)),
    integralOpXView_(integralOpX.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly)),
    alpha_(alpha),
    beta_(beta),
    omega_(omega),
    globalLength_(x.getGlobalLength())
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const std::size_t row) const
  {
    Scalar zero = Kokkos::ArithTraits<Scalar>::zero();
    Scalar one = Kokkos::ArithTraits<Scalar>::one();
    Scalar value = omega_*integralOpXView_(row,0);
    value /= (2*globalLength_*(one-integralOpView_(row,0))*(one-integralOpView_(row,0)));
    if (beta_ == zero)
      yView_(row,0) = alpha_*(value - xView_(row,0));
    else
      yView_(row,0) = alpha_*(value - xView_(row,0)) + beta_*yView_(row,0);
  }

};

//==================================================================

template <class TpetraVectorType, class TpetraMatrixType>
struct PreconditionerEvaluatorFunctor
{
  typedef typename TpetraVectorType::impl_scalar_type Scalar;
  typedef typename TpetraVectorType::local_ordinal_type LO;
  typedef typename TpetraVectorType::global_ordinal_type GO;
  typedef typename TpetraVectorType::execution_space ExecSpace;
  typedef typename TpetraVectorType::dual_view_type::t_dev::const_type ConstViewType;
  typedef typename TpetraMatrixType::local_matrix_device_type LocalMatrixType;

  const LocalMatrixType precLocal_;
  const ConstViewType integralOpView_;
  const Scalar omega_;
  const GO numGlobalElements_;

  PreconditionerEvaluatorFunctor(const TpetraMatrixType& prec,
                                 const TpetraVectorType& integralOp,
                                 const Scalar& omega,
                                 const GO& numGlobalElements) :
    precLocal_(prec.getLocalMatrixDevice()),
    integralOpView_(integralOp.template getLocalView<ExecSpace>(Tpetra::Access::ReadOnly)),
    omega_(omega),
    numGlobalElements_(numGlobalElements)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO localRow) const
  {
    Scalar one = Kokkos::ArithTraits<Scalar>::one();
    Scalar value = one - integralOpView_(localRow, 0);
    value = omega_/(4*numGlobalElements_*value*value) - one;
    value = 1.0/value;
    precLocal_.replaceValues(localRow, &localRow, 1, &value);
  }

};

//==================================================================

template <class TpetraGraphType>
struct GraphSetupFunctor
{
  typedef typename TpetraGraphType::local_graph_device_type::row_map_type::non_const_type OffsetsViewType;
  typedef typename TpetraGraphType::local_graph_device_type::entries_type::non_const_type IndicesViewType;

  const OffsetsViewType offsetsView_;
  const IndicesViewType indicesView_;
  const std::size_t numLocalElements_;

  GraphSetupFunctor(const OffsetsViewType& offsetsView,
                    const IndicesViewType& indicesView,
                    const std::size_t& numLocalElements) :
    offsetsView_(offsetsView),
    indicesView_(indicesView),
    numLocalElements_(numLocalElements)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const std::size_t localRow) const
  {
    offsetsView_(localRow) = localRow;
    if (localRow < numLocalElements_)
      indicesView_(localRow) = localRow;
  }

};

#endif
