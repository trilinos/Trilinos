// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_STENCILPROBLEMS_DECL_HPP
#define GALERI_STENCILPROBLEMS_DECL_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_Problem.hpp"

namespace Galeri::Xpetra {

// =============================================  Scalar Problem =========================================
template <typename Map, typename Matrix, typename MultiVector>
class ScalarProblem : public Problem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;

  ScalarProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<MultiVector> BuildNullspace();
};

// =============================================  Laplace1D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Laplace1DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;

  Laplace1DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
  Teuchos::RCP<RealValuedMultiVector> BuildCoords();
};

// =============================================  Laplace2D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Laplace2DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;

  Laplace2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
  Teuchos::RCP<RealValuedMultiVector> BuildCoords();
};

// =============================================  AnisotropicDiffusion2DProblem  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class AnisotropicDiffusion2DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  AnisotropicDiffusion2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
};

// =============================================  Laplace3D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Laplace3DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;

  Laplace3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
  Teuchos::RCP<RealValuedMultiVector> BuildCoords();
};

// =============================================  Star2D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Star2DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  Star2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
};

// =============================================  BigStar2D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class BigStar2DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  BigStar2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
};

// =============================================  Brick3D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Brick3DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;

  Brick3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
  Teuchos::RCP<RealValuedMultiVector> BuildCoords();
};

// =============================================  Identity  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class IdentityProblem : public Problem<Map, Matrix, MultiVector> {
 public:
  using RealValuedMultiVector = typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector;
  IdentityProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
  Teuchos::RCP<MultiVector> BuildNullspace();
  Teuchos::RCP<RealValuedMultiVector> BuildCoords();
};

// =============================================  Recirc2D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class Recirc2DProblem : public ScalarProblem<Map, Matrix, MultiVector> {
 public:
  Recirc2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map);
  Teuchos::RCP<Matrix> BuildMatrix();
};

}  // namespace Galeri::Xpetra

#endif  // GALERI_STENCILPROBLEMS_DECL_HPP
