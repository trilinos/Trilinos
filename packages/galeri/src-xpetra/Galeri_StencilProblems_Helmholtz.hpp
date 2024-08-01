// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_STENCILPROBLEMS_HELMHOLTZ_HPP
#define GALERI_STENCILPROBLEMS_HELMHOLTZ_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_Problem_Helmholtz.hpp"
#include "Galeri_XpetraMatrixTypes_Helmholtz.hpp"

namespace Galeri {

  namespace Xpetra {

    // =============================================  Helmholtz1D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Helmholtz1DProblem : public Problem_Helmholtz<Map,Matrix,MultiVector> {
    public:
      Helmholtz1DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem_Helmholtz<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix>                                    BuildMatrix();
      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > BuildMatrices();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Helmholtz1DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      double h = this->list_.get("h", 1.0);
      double omega = this->list_.get("omega", 2.0*M_PI);
      double shift = this->list_.get("shift", 0.5);
      Scalar cpxshift(1.0,shift);

      if (nx == -1)
        nx = this->Map_->getGlobalNumElements();

      this->A_ = TriDiag_Helmholtz<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, h, omega, cpxshift);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > Helmholtz1DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrices() {
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      double h = this->list_.get("h", 1.0);
      double omega = this->list_.get("omega", 2.0*M_PI);
      double shift = this->list_.get("shift", 0.5);
      Scalar cpxshift(1.0,shift);

      if (nx == -1)
        nx = this->Map_->getGlobalNumElements();

      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > system;

      system = TriDiag_Helmholtz_Pair<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, h, omega, cpxshift);
      return system;
    }

    // =============================================  Helmholtz2D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Helmholtz2DProblem : public Problem_Helmholtz<Map,Matrix,MultiVector> {
    public:
      Helmholtz2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem_Helmholtz<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix>                                    BuildMatrix();
      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > BuildMatrices();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Helmholtz2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      int model = this->list_.get("model", 0);
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      double h = this->list_.get("h", 1.0);
      double delta = this->list_.get("delta", 1.0);
      int PMLx_left  = this->list_.get("PMLx_left", 0);
      int PMLx_right = this->list_.get("PMLx_right", 0);
      int PMLy_left  = this->list_.get("PMLy_left", 0);
      int PMLy_right = this->list_.get("PMLy_right", 0);
      double omega   = this->list_.get("omega", 2.0*M_PI);
      double shift = this->list_.get("shift", 0.5);
      Scalar cpxshift(1.0,shift);

      if (nx == -1 || ny == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal)sqrt((double)n);
        ny = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
      }

      this->A_ = Cross2D_Helmholtz<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, h, delta,
										 PMLx_left, PMLx_right,
										 PMLy_left, PMLy_right,
										 omega, cpxshift, model);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > Helmholtz2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrices() {
      int model = this->list_.get("model", 0);
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      double h = this->list_.get("h", 1.0);
      double delta = this->list_.get("delta", 1.0);
      int PMLx_left  = this->list_.get("PMLx_left", 0);
      int PMLx_right = this->list_.get("PMLx_right", 0);
      int PMLy_left  = this->list_.get("PMLy_left", 0);
      int PMLy_right = this->list_.get("PMLy_right", 0);
      double omega   = this->list_.get("omega", 2.0*M_PI);
      double shift = this->list_.get("shift", 0.5);
      Scalar cpxshift(1.0,shift);

      if (nx == -1 || ny == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal)sqrt((double)n);
        ny = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
      }

      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > system;

      system = Cross2D_Helmholtz_Pair<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, h, delta,
										    PMLx_left, PMLx_right,
										    PMLy_left, PMLy_right,
										    omega, cpxshift, model);
      return system;
    }

    // =============================================  Helmholtz3D  =============================================
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    class Helmholtz3DProblem : public Problem_Helmholtz<Map,Matrix,MultiVector> {
    public:
      Helmholtz3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem_Helmholtz<Map,Matrix,MultiVector>(list, map) { }
      Teuchos::RCP<Matrix>                                    BuildMatrix();
      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > BuildMatrices();
    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    Teuchos::RCP<Matrix> Helmholtz3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrix() {
      int model = this->list_.get("model", 0);
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      GlobalOrdinal nz = this->list_.get("nz", (GlobalOrdinal) -1);
      double h = this->list_.get("h", 1.0);
      double delta = this->list_.get("delta", 1.0);
      int PMLx_left  = this->list_.get("PMLx_left", 0);
      int PMLx_right = this->list_.get("PMLx_right", 0);
      int PMLy_left  = this->list_.get("PMLy_left", 0);
      int PMLy_right = this->list_.get("PMLy_right", 0);
      int PMLz_left  = this->list_.get("PMLz_left", 0);
      int PMLz_right = this->list_.get("PMLz_right", 0);
      double omega   = this->list_.get("omega", 2.0*M_PI);
      double shift = this->list_.get("shift", 0.5);
      Scalar cpxshift(1.0,shift);

      if (nx == -1 || ny == -1 || nz == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
        ny = nx; nz = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
      }

      this->A_ = Cross3D_Helmholtz<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, h, delta,
										 PMLx_left, PMLx_right,
										 PMLy_left, PMLy_right,
										 PMLz_left, PMLz_right,
										 omega, cpxshift, model);
      this->A_->setObjectLabel(this->getObjectLabel());
      return this->A_;
    }

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > Helmholtz3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>::BuildMatrices() {
      int model = this->list_.get("model", 0);
      GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal) -1);
      GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal) -1);
      GlobalOrdinal nz = this->list_.get("nz", (GlobalOrdinal) -1);
      double h = this->list_.get("h", 1.0);
      double delta = this->list_.get("delta", 1.0);
      int PMLx_left  = this->list_.get("PMLx_left", 0);
      int PMLx_right = this->list_.get("PMLx_right", 0);
      int PMLy_left  = this->list_.get("PMLy_left", 0);
      int PMLy_right = this->list_.get("PMLy_right", 0);
      int PMLz_left  = this->list_.get("PMLz_left", 0);
      int PMLz_right = this->list_.get("PMLz_right", 0);
      double omega   = this->list_.get("omega", 2.0*M_PI);
      double shift = this->list_.get("shift", 0.5);
      Scalar cpxshift(1.0,shift);

      if (nx == -1 || ny == -1 || nz == -1) {
        GlobalOrdinal n = this->Map_->getGlobalNumElements();
        nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
        ny = nx; nz = nx;
        TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
      }

      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > system;

      system = Cross3D_Helmholtz_Pair<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, h, delta,
										    PMLx_left, PMLx_right,
										    PMLy_left, PMLy_right,
										    PMLz_left, PMLz_right,
										    omega, cpxshift, model);
      return system;
    }

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_STENCILPROBLEMS_HELMHOLTZ_HPP
