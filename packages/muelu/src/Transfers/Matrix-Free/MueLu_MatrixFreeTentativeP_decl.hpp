// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
