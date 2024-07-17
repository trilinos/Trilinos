// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_STENCILPROBLEMS_HPP
#define GALERI_STENCILPROBLEMS_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_Problem.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraUtils.hpp"

namespace Galeri {

  namespace Xpetra {

    // =============================================  Scalar Problem =========================================
    template <typename Map, typename Matrix, typename MultiVector>
    class ScalarProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      using RealValuedMultiVector = typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector;

      ScalarProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<MultiVector> BuildNullspace();
    };

    template <typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<MultiVector> ScalarProblem<Map,Matrix,MultiVector>::BuildNullspace() {
      this->Nullspace_ = MultiVectorTraits<Map,MultiVector>::Build(this->Map_, 1);
      this->Nullspace_->putScalar(Teuchos::ScalarTraits<typename MultiVector::scalar_type>::one());
      return this->Nullspace_;
    }

    // =============================================  Laplace1D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Laplace1DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      using RealValuedMultiVector = typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector;

      Laplace1DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
      Teuchos::RCP<RealValuedMultiVector> BuildCoords();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Laplace1DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);

      if (nx == -1)
        nx = this->Map_->getGlobalNumElements();

      this->A_ = TriDiag<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, 2.0, -1.0, -1.0);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector> Laplace1DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildCoords() {

      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;

      if (list.isParameter("nx")) {
        if (list.isType<int>("nx"))
          nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
        else
          nx = list.get<GlobalOrdinal>("nx");
      }

      if (nx == -1) {
        nx = this->Map_->getGlobalNumElements();
      }

      Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("1D", this->Map_, this->list_);

      return this->Coords_;

    }

    // =============================================  Laplace2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Laplace2DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      Laplace2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Laplace2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;

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

      double  one = 1.0;
      Scalar  stretchx = (Scalar) this->list_.get("stretchx", one);
      Scalar  stretchy = (Scalar) this->list_.get("stretchy", one);

      if (nx == -1 || ny == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal)sqrt((double)n);
        ny = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
      }
      bool keepBCs = this->list_.get("keepBCs", false);

      Scalar east   = (Scalar) -one / (stretchx*stretchx);
      Scalar west   = (Scalar) -one / (stretchx*stretchx);
      Scalar north  = (Scalar) -one / (stretchy*stretchy);
      Scalar south  = (Scalar) -one / (stretchy*stretchy);
      Scalar center = -(east + west + north + south);

      this->A_ = Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, center, west, east, south, north, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Laplace2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class AnisotropicDiffusion2DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      AnisotropicDiffusion2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> AnisotropicDiffusion2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
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
      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;
      using MT = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
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
      Scalar one = 1.0;
      Scalar two = one+one;

      MT one_MT  = 1.0;
      MT two_MT = one_MT+one_MT;

      if (nx == -1 || ny == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal)sqrt((double)n);
        ny = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
      }
      bool keepBCs = this->list_.get("keepBCs", false);

      MT   dtInv  = one_MT / this->list_.get("dt", one_MT);
      Scalar Kxx = (Scalar) this->list_.get("Kxx", one);
      Scalar Kxy = (Scalar) this->list_.get("Kxy", zero);
      Scalar Kyy = (Scalar) this->list_.get("Kyy", one);

      std::string meshType = this->list_.get("meshType", "tri");
      Scalar a, b, c, d, e, z1, z2, z3, z4;

      // stencil
      //  z3  e  z4
      //   b  a  c
      //  z1  d  z2

      if (meshType == "tri") {
        MT massDiag_MT = one_MT/two_MT/(Teuchos::as<MT>((nx+1)*(ny+1)));
        Scalar massDiag = massDiag_MT;
        Scalar massOffDiag = massDiag_MT/6.0;

        c  = -Kxx + Kxy + dtInv * massOffDiag;
        b  = -Kxx + Kxy + dtInv * massOffDiag;
        e  = -Kyy + Kxy + dtInv * massOffDiag;
        d  = -Kyy + Kxy + dtInv * massOffDiag;
        a = Kxx * two + Kyy * two - Kxy*two + dtInv * massDiag;
        z1 = -Kxy       + dtInv * massOffDiag;
        z2 = zero;
        z3 = zero;
        z4 = -Kxy       + dtInv * massOffDiag;
      } else if (meshType == "quad") {
        Scalar mass = one/((Scalar) ((nx+1)*(ny+1))) / (Scalar) 36.0;
        Scalar four_thirds = (Scalar) (4.0/3.0);
        Scalar third       = (Scalar) (1.0/3.0);
        Scalar two_thirds  = (Scalar) (2.0/3.0);
        Scalar half        = (Scalar) (1.0/2.0);
        Scalar sixth       = (Scalar) (1.0/6.0);

        a = four_thirds*Kxx  + four_thirds*Kyy             + dtInv * mass * (Scalar) 16.0;

        c  = -two_thirds*Kxx + third*Kyy                   + dtInv * mass * (Scalar) 4.0;
        b  = -two_thirds*Kxx + third*Kyy                   + dtInv * mass * (Scalar) 4.0;
        e  = third*Kxx       - two_thirds*Kyy              + dtInv * mass * (Scalar) 4.0;
        d  = third*Kxx       - two_thirds*Kyy              + dtInv * mass * (Scalar) 4.0;

        z1 = -sixth*Kxx      - sixth*Kyy       - half*Kxy  + dtInv * mass;
        z2 = -sixth*Kxx      - sixth*Kyy       + half*Kxy  + dtInv * mass;
        z3 = -sixth*Kxx      - sixth*Kyy       + half*Kxy  + dtInv * mass;
        z4 = -sixth*Kxx      - sixth*Kyy       - half*Kxy  + dtInv * mass;
      } else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "You need to specify meshType=\"tri\" or \"quad\".");

      this->A_ = Star2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Laplace3D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Laplace3DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      Laplace3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Laplace3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;
      GlobalOrdinal nz = -1;

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
      double  one = 1.0;
      Scalar  stretchx = (Scalar) this->list_.get("stretchx", one);
      Scalar  stretchy = (Scalar) this->list_.get("stretchy", one);
      Scalar  stretchz = (Scalar) this->list_.get("stretchz", one);

      if (nx == -1 || ny == -1 || nz == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
        ny = nx; nz = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
      }
      bool keepBCs = this->list_.get("keepBCs", false);

      Scalar right  = (Scalar) -one / (stretchx*stretchx);
      Scalar left   = (Scalar) -one / (stretchx*stretchx);
      Scalar front  = (Scalar) -one / (stretchy*stretchy);
      Scalar back   = (Scalar) -one / (stretchy*stretchy);
      Scalar up     = (Scalar) -one / (stretchz*stretchz);
      Scalar down   = (Scalar) -one / (stretchz*stretchz);
      Scalar center = -(right + left + front + back + up + down);

      this->A_ = Cross3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, center, left, right, front, back, down, up, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Star2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Star2DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      Star2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Star2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;

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

      Scalar a  = this->list_.get("a",   8.0);
      Scalar b  = this->list_.get("b",  -1.0);
      Scalar c  = this->list_.get("c",  -1.0);
      Scalar d  = this->list_.get("d",  -1.0);
      Scalar e  = this->list_.get("e",  -1.0);
      Scalar z1 = this->list_.get("z1", -1.0);
      Scalar z2 = this->list_.get("z2", -1.0);
      Scalar z3 = this->list_.get("z3", -1.0);
      Scalar z4 = this->list_.get("z4", -1.0);

      bool keepBCs = this->list_.get("keepBCs", false);

      this->A_ = Star2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  BigStar2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class BigStar2DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      BigStar2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> BigStar2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;

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

      this->A_ = BigStar2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, bb, cc, dd, ee, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Brick3D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Brick3DProblem : public ScalarProblem<Map,Matrix,MultiVector> {
    public:
      Brick3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : ScalarProblem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Brick3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;
      GlobalOrdinal nz = -1;

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
        nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
        ny = nx; nz = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
      }

      bool keepBCs = this->list_.get("keepBCs", false);

      this->A_ = Brick3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Identity  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class IdentityProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      using RealValuedMultiVector = typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector;
      IdentityProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
      Teuchos::RCP<MultiVector> BuildNullspace();
      Teuchos::RCP<RealValuedMultiVector> BuildCoords();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> IdentityProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Scalar a = this->list_.get("a", 1.0);
      this->A_ = Identity<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, a);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<MultiVector> IdentityProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildNullspace() {
      this->Nullspace_ = MultiVectorTraits<Map,MultiVector>::Build(this->Map_, 1);
      this->Nullspace_->putScalar(Teuchos::ScalarTraits<typename MultiVector::scalar_type>::one());
      return this->Nullspace_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<typename Problem<Map,Matrix,MultiVector>::RealValuedMultiVector> IdentityProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildCoords() {

      Teuchos::ParameterList list = this->list_;
      GlobalOrdinal nx = -1;

      if (list.isParameter("nx")) {
        if (list.isType<int>("nx"))
          nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
        else
          nx = list.get<GlobalOrdinal>("nx");
      }

      if (nx == -1) {
        nx = this->Map_->getGlobalNumElements();
      }

      Utils::CreateCartesianCoordinates<typename RealValuedMultiVector::scalar_type, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("1D", this->Map_, this->list_);

      return this->Coords_;

    }

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_STENCILPROBLEMS_HPP
