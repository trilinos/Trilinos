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
#ifndef MUELU_GEOMETRICINTERPOLATIONPFACTORY_KOKKOS_DECL_HPP
#define MUELU_GEOMETRICINTERPOLATIONPFACTORY_KOKKOS_DECL_HPP

// Teuchos includes for dense linear algebra
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include "Xpetra_CrsGraph_fwd.hpp"

#include "MueLu_PFactory.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_IndexManager_kokkos_fwd.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class GeometricInterpolationPFactory_kokkos : public PFactory {
#undef MUELU_GEOMETRICINTERPOLATIONPFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  // Declare useful types
  using real_type                  = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using realvaluedmultivector_type = Xpetra::MultiVector<real_type, LO, GO, Node>;
  using device_type                = typename Node::device_type;
  using execution_space            = typename Node::execution_space;
  using impl_scalar_type           = typename Kokkos::ArithTraits<real_type>::val_type;
  using coord_view_type            = typename Kokkos::View<impl_scalar_type**,
                                                Kokkos::LayoutLeft,
                                                device_type>;

  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  GeometricInterpolationPFactory_kokkos() {}

  //! Destructor.
  virtual ~GeometricInterpolationPFactory_kokkos() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Public functors
  //@{

  struct coarseCoordinatesBuilderFunctor {
    IndexManager_kokkos geoData_;
    coord_view_type fineCoordView_;
    coord_view_type coarseCoordView_;

    coarseCoordinatesBuilderFunctor(RCP<IndexManager_kokkos> geoData,
                                    coord_view_type fineCoordView,
                                    coord_view_type coarseCoordView);

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO nodeIdx) const;

  };  // struct coarseCoordinatesBuilderFunctor

  //@}

  void BuildConstantP(RCP<Matrix>& P, RCP<const CrsGraph>& prolongatorGraph, RCP<Matrix>& A) const;

 private:
  void BuildLinearP(RCP<Matrix>& A, RCP<const CrsGraph>& prolongatorGraph,
                    RCP<realvaluedmultivector_type>& fineCoordinates,
                    RCP<realvaluedmultivector_type>& ghostCoordinates,
                    const int numDimensions, RCP<Matrix>& P) const;
  void ComputeLinearInterpolationStencil(const int numDimensions, const int numInterpolationPoints,
                                         const Array<Array<real_type> > coord,
                                         Array<real_type>& stencil) const;
  void GetInterpolationFunctions(const LO numDimensions,
                                 const Teuchos::SerialDenseVector<LO, real_type> parametricCoordinates,
                                 real_type functions[4][8]) const;

};  // class GeometricInterpolationPFactory_kokkos

}  // namespace MueLu

#define MUELU_GEOMETRICINTERPOLATIONPFACTORY_KOKKOS_SHORT
#endif  // MUELU_GEOMETRICINTERPOLATIONPFACTORY_KOKKOS_DECL_HPP
