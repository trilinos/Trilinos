// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef GALERI_STENCILPROBLEMS_HPP
#define GALERI_STENCILPROBLEMS_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_Problem.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"

namespace Galeri {

  namespace Xpetra {

    // =============================================  Laplace1D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Laplace1DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      Laplace1DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
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

    // =============================================  Laplace2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Laplace2DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      Laplace2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Laplace2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      double one = 1.0;
      Scalar  stretchx = this->list_.get("stretchx", one);
      Scalar  stretchy = this->list_.get("stretchy", one);

      if (nx == -1 || ny == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal)sqrt((double)n);
        ny = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
      }
      bool keepBCs = this->list_.get("keepBCs", false);

      Scalar east   = -one / (stretchx*stretchx);
      Scalar west   = -one / (stretchx*stretchx);
      Scalar north  = -one / (stretchy*stretchy);
      Scalar south  = -one / (stretchy*stretchy);
      Scalar center = -(east + west + north + south);

      this->A_ = Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, center, west, east, south, north, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Laplace3D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Laplace3DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      Laplace3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Laplace3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      GlobalOrdinal nz = this->list_.get("nz", (GlobalOrdinal) -1);
      double one = 1.0;
      Scalar  stretchx = this->list_.get("stretchx", one);
      Scalar  stretchy = this->list_.get("stretchy", one);
      Scalar  stretchz = this->list_.get("stretchz", one);

      if (nx == -1 || ny == -1 || nz == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
        ny = nx; nz = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
      }
      bool keepBCs = this->list_.get("keepBCs", false);

      Scalar right  = -one / (stretchx*stretchx);
      Scalar left   = -one / (stretchx*stretchx);
      Scalar front  = -one / (stretchy*stretchy);
      Scalar back   = -one / (stretchy*stretchy);
      Scalar up     = -one / (stretchz*stretchz);
      Scalar down   = -one / (stretchz*stretchz);
      Scalar center = -(right + left + front + back + up + down);

      this->A_ = Cross3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, center, left, right, front, back, down, up, this->DirichletBC_, keepBCs);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Star2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Star2DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      Star2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Star2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);

      Scalar a = this->list_.get("a", 8.0);
      Scalar b = this->list_.get("b", -1.0);
      Scalar c = this->list_.get("c", -1.0);
      Scalar d = this->list_.get("d", -1.0);
      Scalar e = this->list_.get("e", -1.0);
      Scalar z1 = this->list_.get("z1", -1.0);
      Scalar z2 = this->list_.get("z2", -1.0);
      Scalar z3 = this->list_.get("z3", -1.0);
      Scalar z4 = this->list_.get("z4", -1.0);

      this->A_ = Star2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  BigStar2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class BigStar2DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      BigStar2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> BigStar2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);

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

      this->A_ = BigStar2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, a, b, c, d, e, z1, z2, z3, z4, bb, cc, dd, ee);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Brick3D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Brick3DProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      Brick3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Brick3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      GlobalOrdinal nz = this->list_.get("nz", (GlobalOrdinal) -1);

      if (nx == -1 || ny == -1 || nz == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
        ny = nx; nz = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
      }
      this->A_ = Brick3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    // =============================================  Identity  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class IdentityProblem : public Problem<Map,Matrix,MultiVector> {
    public:
      IdentityProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix> BuildMatrix();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> IdentityProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      Scalar a = this->list_.get("a", 1.0);
      this->A_ = Identity<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, a);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_STENCILPROBLEMS_HPP
