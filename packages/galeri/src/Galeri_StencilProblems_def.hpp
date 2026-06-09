// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_STENCILPROBLEMS_DEF_HPP
#define GALERI_STENCILPROBLEMS_DEF_HPP

#include "Galeri_StencilProblems_decl.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraUtils.hpp"

namespace Galeri::Xpetra {

// =============================================  Scalar Problem =========================================
template <typename Map, typename Matrix, typename MultiVector>
ScalarProblem<Map, Matrix, MultiVector>::ScalarProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : Problem<Map, Matrix, MultiVector>(list, map) {}

template <typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<MultiVector> ScalarProblem<Map, Matrix, MultiVector>::BuildNullspace() {
  this->Nullspace_ = MultiVectorTraits<Map, MultiVector>::Build(this->Map_, 1);
  this->Nullspace_->putScalar(Teuchos::ScalarTraits<typename MultiVector::scalar_type>::one());
  return this->Nullspace_;
}

// =============================================  Laplace1D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Laplace1DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Laplace1DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Laplace1DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal)-1);

  if (nx == -1)
    nx = this->Map_->getGlobalNumElements();

  bool keepBCs = false;
  this->A_     = TriDiag<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, 2.0, -1.0, -1.0, this->DirichletBC_, keepBCs, "Laplace 1D");
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> Laplace1DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }

  if (nx == -1) {
    nx = this->Map_->getGlobalNumElements();
  }

  this->Coords_ = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("1D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  Laplace2D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Laplace2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Laplace2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Laplace2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }

  double one      = 1.0;
  Scalar stretchx = (Scalar)this->list_.get("stretchx", one);
  Scalar stretchy = (Scalar)this->list_.get("stretchy", one);

  if (nx == -1 || ny == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)sqrt((double)n);
    ny              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny != n, std::logic_error, "You need to specify nx and ny.");
  }
  bool keepBCs = this->list_.get("keepBCs", false);

  Scalar east   = (Scalar)-one / (stretchx * stretchx);
  Scalar west   = (Scalar)-one / (stretchx * stretchx);
  Scalar north  = (Scalar)-one / (stretchy * stretchy);
  Scalar south  = (Scalar)-one / (stretchy * stretchy);
  Scalar center = -(east + west + north + south);

  this->A_ = Cross2D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, center, west, east, south, north, this->DirichletBC_, keepBCs, "Laplace 2D");
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> Laplace2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  this->Coords_               = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("2D", this->Map_, this->list_);
  return this->Coords_;
}

// =============================================  AnisotropicDiffusion2DProblem  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
AnisotropicDiffusion2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::AnisotropicDiffusion2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> AnisotropicDiffusion2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  /*
    2D triangle mesh (meshType == tri)

    *-*-*
    |/|/|
    *-*-*
    |/|/|
    *-*-*

    or 2D triangle mesh (meshType == quad)

    *-*-*
    | | |
    *-*-*
    | | |
    *-*-*

    using piecewise linear continuous elements

    PDE

     -\nabla K \nabla u + 1/dt * u = f

    for some constant symmetric 2x2 matrix

    K = [[Kxx, Kxy],
         [Kxy, Kyy]]
  */

  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;
  using MT                    = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }

  Scalar zero = 0.0;
  Scalar one  = 1.0;
  Scalar two  = one + one;

  MT one_MT = 1.0;
  MT two_MT = one_MT + one_MT;

  if (nx == -1 || ny == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)sqrt((double)n);
    ny              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny != n, std::logic_error, "You need to specify nx and ny.");
  }
  bool keepBCs = this->list_.get("keepBCs", false);

  MT dtInv   = one_MT / this->list_.get("dt", one_MT);
  Scalar Kxx = (Scalar)this->list_.get("Kxx", one);
  Scalar Kxy = (Scalar)this->list_.get("Kxy", zero);
  Scalar Kyy = (Scalar)this->list_.get("Kyy", one);

  std::string meshType = this->list_.get("meshType", "tri");
  Scalar a, b, c, d, e, z1, z2, z3, z4;

  // stencil
  //  z3  e  z4
  //   b  a  c
  //  z1  d  z2

  if (meshType == "tri") {
    MT massDiag_MT     = one_MT / two_MT / (Teuchos::as<MT>((nx + 1) * (ny + 1)));
    Scalar massDiag    = massDiag_MT;
    Scalar massOffDiag = massDiag_MT / 6.0;

    c  = -Kxx + Kxy + dtInv * massOffDiag;
    b  = -Kxx + Kxy + dtInv * massOffDiag;
    e  = -Kyy + Kxy + dtInv * massOffDiag;
    d  = -Kyy + Kxy + dtInv * massOffDiag;
    a  = Kxx * two + Kyy * two - Kxy * two + dtInv * massDiag;
    z1 = -Kxy + dtInv * massOffDiag;
    z2 = zero;
    z3 = zero;
    z4 = -Kxy + dtInv * massOffDiag;
  } else if (meshType == "quad") {
    Scalar mass        = one / ((Scalar)((nx + 1) * (ny + 1))) / (Scalar)36.0;
    Scalar four_thirds = (Scalar)(4.0 / 3.0);
    Scalar third       = (Scalar)(1.0 / 3.0);
    Scalar two_thirds  = (Scalar)(2.0 / 3.0);
    Scalar half        = (Scalar)(1.0 / 2.0);
    Scalar sixth       = (Scalar)(1.0 / 6.0);

    a = four_thirds * Kxx + four_thirds * Kyy + dtInv * mass * (Scalar)16.0;

    c = -two_thirds * Kxx + third * Kyy + dtInv * mass * (Scalar)4.0;
    b = -two_thirds * Kxx + third * Kyy + dtInv * mass * (Scalar)4.0;
    e = third * Kxx - two_thirds * Kyy + dtInv * mass * (Scalar)4.0;
    d = third * Kxx - two_thirds * Kyy + dtInv * mass * (Scalar)4.0;

    z1 = -sixth * Kxx - sixth * Kyy - half * Kxy + dtInv * mass;
    z2 = -sixth * Kxx - sixth * Kyy + half * Kxy + dtInv * mass;
    z3 = -sixth * Kxx - sixth * Kyy + half * Kxy + dtInv * mass;
    z4 = -sixth * Kxx - sixth * Kyy - half * Kxy + dtInv * mass;
  } else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "You need to specify meshType=\"tri\" or \"quad\".");

  this->A_ = Star2D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, this->DirichletBC_, keepBCs);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

// =============================================  Laplace3D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Laplace3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Laplace3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Laplace3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;
  GlobalOrdinal nz            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }
  if (list.isParameter("nz")) {
    if (list.isType<int>("nz"))
      nz = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
    else
      nz = list.get<GlobalOrdinal>("nz");
  }
  double one      = 1.0;
  Scalar stretchx = (Scalar)this->list_.get("stretchx", one);
  Scalar stretchy = (Scalar)this->list_.get("stretchy", one);
  Scalar stretchz = (Scalar)this->list_.get("stretchz", one);

  if (nx == -1 || ny == -1 || nz == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny              = nx;
    nz              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }
  bool keepBCs = this->list_.get("keepBCs", false);

  Scalar right  = (Scalar)-one / (stretchx * stretchx);
  Scalar left   = (Scalar)-one / (stretchx * stretchx);
  Scalar front  = (Scalar)-one / (stretchy * stretchy);
  Scalar back   = (Scalar)-one / (stretchy * stretchy);
  Scalar up     = (Scalar)-one / (stretchz * stretchz);
  Scalar down   = (Scalar)-one / (stretchz * stretchz);
  Scalar center = -(right + left + front + back + up + down);

  this->A_ = Cross3D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, nz, center, left, right, front, back, down, up, this->DirichletBC_, keepBCs, "Laplace 3D");
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> Laplace3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  this->Coords_               = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("3D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  Star2D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Star2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Star2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Star2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }

  Scalar a  = this->list_.get("a", 8.0);
  Scalar b  = this->list_.get("b", -1.0);
  Scalar c  = this->list_.get("c", -1.0);
  Scalar d  = this->list_.get("d", -1.0);
  Scalar e  = this->list_.get("e", -1.0);
  Scalar z1 = this->list_.get("z1", -1.0);
  Scalar z2 = this->list_.get("z2", -1.0);
  Scalar z3 = this->list_.get("z3", -1.0);
  Scalar z4 = this->list_.get("z4", -1.0);

  bool keepBCs = this->list_.get("keepBCs", false);

  this->A_ = Star2D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, this->DirichletBC_, keepBCs);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

// =============================================  BigStar2D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
BigStar2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BigStar2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> BigStar2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }

  Scalar a  = this->list_.get("a", 20.0);
  Scalar b  = this->list_.get("b", -8.0);
  Scalar c  = this->list_.get("c", -8.0);
  Scalar d  = this->list_.get("d", -8.0);
  Scalar e  = this->list_.get("e", -8.0);
  Scalar z1 = this->list_.get("z1", 2.0);
  Scalar z2 = this->list_.get("z2", 2.0);
  Scalar z3 = this->list_.get("z3", 2.0);
  Scalar z4 = this->list_.get("z4", 2.0);
  Scalar bb = this->list_.get("bb", 1.0);
  Scalar cc = this->list_.get("cc", 1.0);
  Scalar dd = this->list_.get("dd", 1.0);
  Scalar ee = this->list_.get("ee", 1.0);

  bool keepBCs = this->list_.get("keepBCs", false);

  this->A_ = BigStar2D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, bb, cc, dd, ee, this->DirichletBC_, keepBCs);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

// =============================================  Brick3D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Brick3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Brick3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Brick3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;
  GlobalOrdinal nz            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }
  if (list.isParameter("nz")) {
    if (list.isType<int>("nz"))
      nz = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
    else
      nz = list.get<GlobalOrdinal>("nz");
  }

  if (nx == -1 || ny == -1 || nz == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny              = nx;
    nz              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }

  bool keepBCs = this->list_.get("keepBCs", false);

  this->A_ = Brick3D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, this->DirichletBC_, keepBCs, "3D 27 point stencil");
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> Brick3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  this->Coords_               = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("3D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  Scalar3D_27Pt =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Scalar3D_27PtProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Scalar3D_27PtProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Scalar3D_27PtProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;
  GlobalOrdinal nz            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }
  if (list.isParameter("nz")) {
    if (list.isType<int>("nz"))
      nz = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
    else
      nz = list.get<GlobalOrdinal>("nz");
  }

  if (nx == -1 || ny == -1 || nz == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny              = nx;
    nz              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }
  Scalar S111 = -1.0, S211 = -1.0, S311 = -1.0;
  Scalar S121 = -1.0, S221 = -1.0, S321 = -1.0;
  Scalar S131 = -1.0, S231 = -1.0, S331 = -1.0;
  Scalar S112 = -1.0, S212 = -1.0, S312 = -1.0;
  Scalar S122 = -1.0, S222 = 26.0, S322 = -1.0;
  Scalar S132 = -1.0, S232 = -1.0, S332 = -1.0;
  Scalar S113 = -1.0, S213 = -1.0, S313 = -1.0;
  Scalar S123 = -1.0, S223 = -1.0, S323 = -1.0;
  Scalar S133 = -1.0, S233 = -1.0, S333 = -1.0;

  // 27-pt stencil given by
  //
  //            ^    S131  S231 S331
  //     z=-1   |y   S121  S221 S321
  //            |    S111  S211 S311
  //                 ----- x ------>
  //
  //            ^    S132  S232 S332
  //     z= 0   |y   S122  S222 S322
  //            |    S112  S212 S312
  //                 ----- x ------>
  //
  //            ^    S133  S233 S333
  //     z=+1   |y   S123  S223 S323
  //            |    S113  S213 S313
  //                 ----- x ------>
  //
  if ((list.isParameter("S111")) && (list.isType<Scalar>("S111"))) S111 = list.get("S111", -1);
  if ((list.isParameter("S112")) && (list.isType<Scalar>("S112"))) S112 = list.get("S112", -1);
  if ((list.isParameter("S113")) && (list.isType<Scalar>("S113"))) S113 = list.get("S113", -1);
  if ((list.isParameter("S121")) && (list.isType<Scalar>("S121"))) S121 = list.get("S121", -1);
  if ((list.isParameter("S122")) && (list.isType<Scalar>("S122"))) S122 = list.get("S122", -1);
  if ((list.isParameter("S123")) && (list.isType<Scalar>("S123"))) S123 = list.get("S123", -1);
  if ((list.isParameter("S131")) && (list.isType<Scalar>("S131"))) S131 = list.get("S131", -1);
  if ((list.isParameter("S132")) && (list.isType<Scalar>("S132"))) S132 = list.get("S132", -1);
  if ((list.isParameter("S133")) && (list.isType<Scalar>("S133"))) S133 = list.get("S133", -1);

  if ((list.isParameter("S211")) && (list.isType<Scalar>("S211"))) S211 = list.get("S211", -1);
  if ((list.isParameter("S212")) && (list.isType<Scalar>("S212"))) S212 = list.get("S212", -1);
  if ((list.isParameter("S213")) && (list.isType<Scalar>("S213"))) S213 = list.get("S213", -1);
  if ((list.isParameter("S221")) && (list.isType<Scalar>("S221"))) S221 = list.get("S221", -1);
  if ((list.isParameter("S222")) && (list.isType<Scalar>("S222"))) S222 = list.get("S222", 26.0);
  if ((list.isParameter("S223")) && (list.isType<Scalar>("S223"))) S223 = list.get("S223", -1);
  if ((list.isParameter("S231")) && (list.isType<Scalar>("S231"))) S231 = list.get("S231", -1);
  if ((list.isParameter("S232")) && (list.isType<Scalar>("S232"))) S232 = list.get("S232", -1);
  if ((list.isParameter("S233")) && (list.isType<Scalar>("S233"))) S233 = list.get("S233", -1);

  if ((list.isParameter("S311")) && (list.isType<Scalar>("S311"))) S311 = list.get("S311", -1);
  if ((list.isParameter("S312")) && (list.isType<Scalar>("S312"))) S312 = list.get("S312", -1);
  if ((list.isParameter("S313")) && (list.isType<Scalar>("S313"))) S313 = list.get("S313", -1);
  if ((list.isParameter("S321")) && (list.isType<Scalar>("S321"))) S321 = list.get("S321", -1);
  if ((list.isParameter("S322")) && (list.isType<Scalar>("S322"))) S322 = list.get("S322", -1);
  if ((list.isParameter("S323")) && (list.isType<Scalar>("S323"))) S323 = list.get("S323", -1);
  if ((list.isParameter("S331")) && (list.isType<Scalar>("S331"))) S331 = list.get("S331", -1);
  if ((list.isParameter("S332")) && (list.isType<Scalar>("S332"))) S332 = list.get("S332", -1);
  if ((list.isParameter("S333")) && (list.isType<Scalar>("S333"))) S333 = list.get("S333", -1);

  bool keepBCs = this->list_.get("keepBCs", false);

  this->A_ = Scalar3D_27Pt<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, nz,
                                                                             S111, S211, S311, S121, S221, S321, S131, S231, S331,
                                                                             S112, S212, S312, S122, S222, S322, S132, S232, S332,
                                                                             S113, S213, S313, S123, S223, S323, S133, S233, S333,
                                                                             this->DirichletBC_, keepBCs, "3D 27 point stencil", false);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> Scalar3D_27PtProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  this->Coords_               = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("3D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  HexFEM_LapStiff =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
HexFEM_LapStiffProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::HexFEM_LapStiffProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> HexFEM_LapStiffProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;
  GlobalOrdinal nz            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }
  if (list.isParameter("nz")) {
    if (list.isType<int>("nz"))
      nz = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
    else
      nz = list.get<GlobalOrdinal>("nz");
  }

  if (nx == -1 || ny == -1 || nz == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny              = nx;
    nz              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }
  Scalar one      = Teuchos::ScalarTraits<Scalar>::one();
  Scalar stretchx = (Scalar)this->list_.get("stretchx", one);
  Scalar stretchy = (Scalar)this->list_.get("stretchy", one);
  Scalar stretchz = (Scalar)this->list_.get("stretchz", one);
  Scalar lx       = (Scalar)this->list_.get("lx", one);
  Scalar ly       = (Scalar)this->list_.get("ly", one);
  Scalar lz       = (Scalar)this->list_.get("lz", one);
  Scalar hx, hy, hz;

  // not sure why in other galeri spots nx+1,ny+1,nz+1 are used. When
  // coordinates printed, nx-1,ny-1,nz-1 must be used to match coordinates
  hx = stretchx * lx / ((Scalar)(nx - 1));
  hy = stretchy * ly / ((Scalar)(ny - 1));
  hz = stretchz * lz / ((Scalar)(nz - 1));

  // Unassembled stiffness matrix coefficients

  Scalar Sdiag     = (hy * hz) / (9. * hx) + (hx * (hy * hy + hz * hz)) / (9. * hy * hz);
  Scalar Sxneigh   = (-1. / 9.) * (hy * hz) / hx + (hx * (hy * hy + hz * hz)) / (18. * hy * hz);
  Scalar Syneigh   = (hx * hy) / (18. * hz) - (hx * hz) / (9. * hy) + (hy * hz) / (18. * hx);
  Scalar Szneigh   = (-1 / 9.) * (hx * hy) / hz + (hx * hz) / (18. * hy) + (hy * hz) / (18. * hx);
  Scalar Sxyneigh  = (hx * hy) / (36. * hz) - (hx * hz) / (18. * hy) - (hy * hz) / (18. * hx);
  Scalar Sxzneigh  = (-1 / 18.) * (hx * hy) / hz + (hx * hz) / (36. * hy) - (hy * hz) / (18. * hx);
  Scalar Syzneigh  = (hy * hz) / (36. * hx) - (hx * (hy * hy + hz * hz)) / (18. * hy * hz);
  Scalar Sxyzneigh = (-1 / 36.) * (hy * hz) / hx - (hx * (hy * hy + hz * hz)) / (36. * hy * hz);

  // We assume an interior stencil. So the diagonal entry has 8
  // contributions, while off-diagonals have either 1, 2, or 4
  // contributions
  Scalar S111 = Sxyzneigh, S211 = Syzneigh * 2., S311 = Sxyzneigh;
  Scalar S121 = Sxzneigh * 2.0, S221 = Szneigh * 4., S321 = Sxzneigh * 2.0;
  Scalar S131 = Sxyzneigh, S231 = Syzneigh * 2., S331 = Sxyzneigh;

  Scalar S112 = Sxyneigh * 2., S212 = Syneigh * 4., S312 = Sxyneigh * 2.;
  Scalar S122 = Sxneigh * 4., S222 = Sdiag * 8., S322 = Sxneigh * 4.;
  Scalar S132 = Sxyneigh * 2., S232 = Syneigh * 4., S332 = Sxyneigh * 2.;

  Scalar S113 = Sxyzneigh, S213 = Syzneigh * 2., S313 = Sxyzneigh;
  Scalar S123 = Sxzneigh * 2.0, S223 = Szneigh * 4., S323 = Sxzneigh * 2.0;
  Scalar S133 = Sxyzneigh, S233 = Syzneigh * 2., S333 = Sxyzneigh;

  // 27-pt stencil given by
  //
  //            ^    S131  S231 S331
  //     z=-1   |y   S121  S221 S321
  //            |    S111  S211 S311
  //                 ----- x ------>
  //
  //            ^    S132  S232 S332
  //     z= 0   |y   S122  S222 S322
  //            |    S112  S212 S312
  //                 ----- x ------>
  //
  //            ^    S133  S233 S333
  //     z=+1   |y   S123  S223 S323
  //            |    S113  S213 S313
  //                 ----- x ------>
  //

  bool keepBCs = this->list_.get("keepBCs", false);

  this->A_ = Scalar3D_27Pt<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, nz,
                                                                             S111, S211, S311, S121, S221, S321, S131, S231, S331,
                                                                             S112, S212, S312, S122, S222, S322, S132, S232, S332,
                                                                             S113, S213, S313, S123, S223, S323, S133, S233, S333,
                                                                             this->DirichletBC_, keepBCs, "HexFEM_LapStiff", true);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> HexFEM_LapStiffProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  this->Coords_               = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("3D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  HexFEM_Mass =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
HexFEM_MassProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::HexFEM_MassProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> HexFEM_MassProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;
  GlobalOrdinal nz            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }
  if (list.isParameter("nz")) {
    if (list.isType<int>("nz"))
      nz = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
    else
      nz = list.get<GlobalOrdinal>("nz");
  }

  if (nx == -1 || ny == -1 || nz == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny              = nx;
    nz              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }
  // compute mesh spacing
  Scalar one      = Teuchos::ScalarTraits<Scalar>::one();
  Scalar stretchx = (Scalar)this->list_.get("stretchx", one);
  Scalar stretchy = (Scalar)this->list_.get("stretchy", one);
  Scalar stretchz = (Scalar)this->list_.get("stretchz", one);
  Scalar lx       = (Scalar)this->list_.get("lx", one);
  Scalar ly       = (Scalar)this->list_.get("ly", one);
  Scalar lz       = (Scalar)this->list_.get("lz", one);
  Scalar hx, hy, hz;

  // not sure why in other galeri spots nx+1,ny+1,nz+1 are used. When
  // coordinates printed, nx-1,ny-1,nz-1 must be used to match coordinates
  hx = stretchx * lx / ((Scalar)(nx - 1));
  hy = stretchy * ly / ((Scalar)(ny - 1));
  hz = stretchz * lz / ((Scalar)(nz - 1));

  Scalar alpha = hx * hy * hz / 216.;

  Scalar M111 = alpha * 1.0, M211 = alpha * 4.0, M311 = alpha * 1.0;
  Scalar M121 = alpha * 4.0, M221 = alpha * 16.0, M321 = alpha * 4.0;
  Scalar M131 = alpha * 1.0, M231 = alpha * 4.0, M331 = alpha * 1.0;

  Scalar M112 = alpha * 4.0, M212 = alpha * 16.0, M312 = alpha * 4.0;
  Scalar M122 = alpha * 16.0, M222 = alpha * 64.0, M322 = alpha * 16.0;
  Scalar M132 = alpha * 4.0, M232 = alpha * 16.0, M332 = alpha * 4.0;

  Scalar M113 = alpha * 1.0, M213 = alpha * 4.0, M313 = alpha * 1.0;
  Scalar M123 = alpha * 4.0, M223 = alpha * 16.0, M323 = alpha * 4.0;
  Scalar M133 = alpha * 1.0, M233 = alpha * 4.0, M333 = alpha * 1.0;

  // 27-pt stencil given by
  //
  //            ^    M131  M231 M331
  //     z=-1   |y   M121  M221 M321
  //            |    M111  M211 M311
  //                 ----- x ------>
  //
  //            ^    M132  M232 M332
  //     z= 0   |y   M122  M222 M322
  //            |    M112  M212 M312
  //                 ----- x ------>
  //
  //            ^    M133  M233 M333
  //     z=+1   |y   M123  M223 M323
  //            |    M113  M213 M313
  //                 ----- x ------>
  //

  bool keepBCs = this->list_.get("keepBCs", false);

  this->A_ = Scalar3D_27Pt<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, nz,
                                                                             M111, M211, M311, M121, M221, M321, M131, M231, M331,
                                                                             M112, M212, M312, M122, M222, M322, M132, M232, M332,
                                                                             M113, M213, M313, M123, M223, M323, M133, M233, M333,
                                                                             this->DirichletBC_, keepBCs, "HexFEM_Mass", true);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> HexFEM_MassProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  this->Coords_               = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("3D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  Identity  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
IdentityProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::IdentityProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : Problem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> IdentityProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Scalar a = this->list_.get("a", 1.0);

  GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal)-1);

  if (nx == -1)
    nx = this->Map_->getGlobalNumElements();

  this->A_ = ScaledIdentity<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, a, "Scaled Identity");
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<MultiVector> IdentityProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildNullspace() {
  this->Nullspace_ = MultiVectorTraits<Map, MultiVector>::Build(this->Map_, 1);
  this->Nullspace_->putScalar(Teuchos::ScalarTraits<typename MultiVector::scalar_type>::one());
  return this->Nullspace_;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<typename Problem<Map, Matrix, MultiVector>::RealValuedMultiVector> IdentityProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildCoords() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;

  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }

  if (nx == -1) {
    nx = this->Map_->getGlobalNumElements();
  }

  this->Coords_ = Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("1D", this->Map_, this->list_);

  return this->Coords_;
}

// =============================================  Recirc2D  =============================================

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Recirc2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::Recirc2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
  : ScalarProblem<Map, Matrix, MultiVector>(list, map) {}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> Recirc2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  Teuchos::ParameterList list = this->list_;
  GlobalOrdinal nx            = -1;
  GlobalOrdinal ny            = -1;

  // TODO FIXME need to check that all parameters are provided, otherwise error out
  if (list.isParameter("nx")) {
    if (list.isType<int>("nx"))
      nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
    else
      nx = list.get<GlobalOrdinal>("nx");
  }
  if (list.isParameter("ny")) {
    if (list.isType<int>("ny"))
      ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
    else
      ny = list.get<GlobalOrdinal>("ny");
  }

  double lx = 1.;
  if (list.isParameter("lx")) {
    lx = list.get<double>("lx");
  }
  double ly = 1.;
  if (list.isParameter("ly")) {
    ly = list.get<double>("ly");
  }

  double conv = 1.;
  if (list.isParameter("convection")) {
    conv = list.get<double>("convection");
  }

  double diff = 1.;
  if (list.isParameter("diffusion")) {
    diff = list.get<double>("diffusion");
  }

  auto& map                                                = this->Map_;
  LocalOrdinal numMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  auto A = MultiVectorTraits<Map, MultiVector>::Build(map, 1);
  auto B = MultiVectorTraits<Map, MultiVector>::Build(map, 1);
  auto C = MultiVectorTraits<Map, MultiVector>::Build(map, 1);
  auto D = MultiVectorTraits<Map, MultiVector>::Build(map, 1);
  auto E = MultiVectorTraits<Map, MultiVector>::Build(map, 1);

  double zero = Teuchos::ScalarTraits<double>::zero();

  A->putScalar(zero);
  B->putScalar(zero);
  C->putScalar(zero);
  D->putScalar(zero);
  E->putScalar(zero);

  auto Adata = A->getDataNonConst(0);
  auto Bdata = B->getDataNonConst(0);
  auto Cdata = C->getDataNonConst(0);
  auto Ddata = D->getDataNonConst(0);
  auto Edata = E->getDataNonConst(0);

  double hx = lx / (nx + 1);
  double hy = ly / (ny + 1);

  for (int i = 0; i < numMyElements; i++) {
    int ix, iy;
    ix           = (myGlobalElements[i]) % nx;
    iy           = (myGlobalElements[i] - ix) / nx;
    double x     = hx * (ix + 1);
    double y     = hy * (iy + 1);
    double ConvX = conv * 4 * x * (x - 1.) * (1. - 2 * y) / hx;
    double ConvY = -conv * 4 * y * (y - 1.) * (1. - 2 * x) / hy;

    // convection part

    if (ConvX < zero) {
      Cdata[i] += ConvX;
      Adata[i] -= ConvX;
    } else {
      Bdata[i] -= ConvX;
      Adata[i] += ConvX;
    }

    if (ConvY < zero) {
      Edata[i] += ConvY;
      Adata[i] -= ConvY;
    } else {
      Ddata[i] -= ConvY;
      Adata[i] += ConvY;
    }

    // add diffusion part
    Adata[i] += diff * 2. / (hx * hx) + diff * 2. / (hy * hy);
    Bdata[i] -= diff / (hx * hx);
    Cdata[i] -= diff / (hx * hx);
    Ddata[i] -= diff / (hy * hy);
    Edata[i] -= diff / (hy * hy);
  }
  Adata = Bdata = Cdata = Ddata = Edata = Teuchos::null;

  this->A_ = Cross2D<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix>(this->Map_, nx, ny, A, B, C, D, E);
  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}  // Recirc2DProblem

}  // namespace Galeri::Xpetra

#define GALERI_STENCILPROBLEMS_INSTANT_TPETRA(S, LO, GO, N)                                                                                                            \
  template class Galeri::Xpetra::ScalarProblem<Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;                            \
  template class Galeri::Xpetra::Laplace1DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::Laplace2DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::AnisotropicDiffusion2DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>; \
  template class Galeri::Xpetra::Laplace3DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::Star2DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;                 \
  template class Galeri::Xpetra::BigStar2DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::Brick3DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;                \
  template class Galeri::Xpetra::Scalar3D_27PtProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;          \
  template class Galeri::Xpetra::HexFEM_LapStiffProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;        \
  template class Galeri::Xpetra::HexFEM_MassProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;            \
  template class Galeri::Xpetra::IdentityProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;               \
  template class Galeri::Xpetra::Recirc2DProblem<S, LO, GO, Tpetra::Map<LO, GO, N>, Tpetra::CrsMatrix<S, LO, GO, N>, Tpetra::MultiVector<S, LO, GO, N>>;

#define GALERI_STENCILPROBLEMS_INSTANT_XPETRA(S, LO, GO, N)                                                                                                                  \
  template class Galeri::Xpetra::ScalarProblem<Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                              \
  template class Galeri::Xpetra::Laplace1DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                \
  template class Galeri::Xpetra::Laplace2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                \
  template class Galeri::Xpetra::AnisotropicDiffusion2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;   \
  template class Galeri::Xpetra::Laplace3DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                \
  template class Galeri::Xpetra::Star2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                   \
  template class Galeri::Xpetra::BigStar2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                \
  template class Galeri::Xpetra::Brick3DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                  \
  template class Galeri::Xpetra::Scalar3D_27PtProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;            \
  template class Galeri::Xpetra::HexFEM_LapStiffProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;          \
  template class Galeri::Xpetra::HexFEM_MassProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::IdentityProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                 \
  template class Galeri::Xpetra::Recirc2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::CrsMatrixWrap<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                 \
                                                                                                                                                                             \
  template class Galeri::Xpetra::ScalarProblem<Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                            \
  template class Galeri::Xpetra::Laplace1DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::Laplace2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::AnisotropicDiffusion2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>; \
  template class Galeri::Xpetra::Laplace3DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::Star2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                 \
  template class Galeri::Xpetra::BigStar2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;              \
  template class Galeri::Xpetra::Brick3DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;                \
  template class Galeri::Xpetra::Scalar3D_27PtProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;          \
  template class Galeri::Xpetra::HexFEM_LapStiffProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;        \
  template class Galeri::Xpetra::HexFEM_MassProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;            \
  template class Galeri::Xpetra::IdentityProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;               \
  template class Galeri::Xpetra::Recirc2DProblem<S, LO, GO, Xpetra::Map<LO, GO, N>, Xpetra::TpetraCrsMatrix<S, LO, GO, N>, Xpetra::MultiVector<S, LO, GO, N>>;

#ifdef HAVE_GALERI_XPETRA

#define GALERI_STENCILPROBLEMS_INSTANT(S, LO, GO, N)  \
  GALERI_STENCILPROBLEMS_INSTANT_TPETRA(S, LO, GO, N) \
  GALERI_STENCILPROBLEMS_INSTANT_XPETRA(S, LO, GO, N)

#else

#define GALERI_STENCILPROBLEMS_INSTANT(S, LO, GO, N) \
  GALERI_STENCILPROBLEMS_INSTANT_TPETRA(S, LO, GO, N)

#endif

#endif
