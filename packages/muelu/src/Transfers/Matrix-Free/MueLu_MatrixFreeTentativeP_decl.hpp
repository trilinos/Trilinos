// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATRIXFREETENTATIVEP_DECL_HPP
#define MUELU_MATRIXFREETENTATIVEP_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_MatrixFreeTentativeP_fwd.hpp"

#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Teuchos_BLAS_types.hpp>
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_Operator_fwd.hpp"

namespace MueLu {

/*!
  @class MatrixFreeTentativeP class.
  @brief Matrix-free tentative restrictor operator.
*/
// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
// class MatrixFreeTentativeP;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MatrixFreeTentativeP : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef typename Node::execution_space execution_space;
  typedef Kokkos::RangePolicy<local_ordinal_type, execution_space> range_type;
  typedef Kokkos::MDRangePolicy<local_ordinal_type, execution_space, Kokkos::Rank<2>> md_range_type;
  typedef Node node_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType real_type;

 private:
#undef MUELU_MATRIXFREETENTATIVEP_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  MatrixFreeTentativeP(Teuchos::RCP<const Map> coarse_map, Teuchos::RCP<const Map> fine_map, Teuchos::RCP<const Aggregates> aggregates)
    : fine_map_(fine_map)
    , coarse_map_(coarse_map)
    , aggregates_(aggregates) {}

  //! Destructor.
  ~MatrixFreeTentativeP() = default;
  //@}

  // compute the apply operator, Y = alpha*R*X + beta*Y
  void apply(const MultiVector &X, MultiVector &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(), Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override;

  // compute the residual
  void residual(const MultiVector &X, const MultiVector &B, MultiVector &R) const override;

  // get the range map
  const Teuchos::RCP<const Map> getRangeMap() const override {
    return fine_map_;
  }

  // get the domain map
  const Teuchos::RCP<const Map> getDomainMap() const override {
    return coarse_map_;
  }

  // get the aggregates
  Teuchos::RCP<const Aggregates> getAggregates() const {
    return aggregates_;
  }

 private:
  // the fine map
  const Teuchos::RCP<const Map> fine_map_;

  // the coarse map
  const Teuchos::RCP<const Map> coarse_map_;

  // the aggregates required for the grid transfer
  const Teuchos::RCP<const Aggregates> aggregates_;
};

}  // namespace MueLu

#define MUELU_MATRIXFREETENTATIVEP_SHORT
#endif  // MUELU_MATRIXFREETENTATIVEP_DECL_HPP
