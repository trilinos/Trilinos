// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef CREATEADRMATRIX_HPP
#define CREATEADRMATRIX_HPP

#include <Teuchos_RCP.hpp>

#include "ADR_Problem.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"

namespace ADR {

namespace Xpetra {

// =============================================  ADR1D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class ADR1DProblem : public Problem<Map, Matrix, MultiVector> {
 public:
  ADR1DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
    : Problem<Map, Matrix, MultiVector>(list, map) {}
  Teuchos::RCP<Matrix> BuildMatrix();

 private:
  // domain definition
  Scalar xleft  = 0.0;
  Scalar xright = 1.0;

  // Function that defines the diffusion coefficient
  inline Scalar diff(Scalar x) { return 1.0 + x; };

  // Function that defines the first derivative of the diffusion coefficient
  inline Scalar diff_prime(Scalar x) { return 1.0; };

  // Function that defines the advection coefficient
  inline Scalar adv(Scalar x) { return 10.0 * x; };

  // Function taht defines the reaction coefficient
  inline Scalar reac(Scalar x) { return 0.0; };
};

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> ADR1DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal)-1);

  if (nx == -1)
    nx = this->Map_->getGlobalNumElements();

  const Scalar dx = (Scalar)(xright - xleft) / static_cast<Scalar>(nx - 1);

  const Scalar a = 2.0;
  const Scalar b = -1.0;
  const Scalar c = -1.0;

  // this->A_ = Galeri::Xpetra::TriDiag<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, 2.0, -1.0, -1.0);
  Teuchos::RCP<Matrix> mtx = Galeri::Xpetra::MatrixTraits<Map, Matrix>::Build(this->Map_, 3);

  LocalOrdinal NumMyElements                               = this->Map_->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = this->Map_->getLocalElementList();
  GlobalOrdinal indexBase                                  = this->Map_->getIndexBase();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = this->Map_->getComm();

  GlobalOrdinal NumGlobalElements = this->Map_->getGlobalNumElements();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  comm->barrier();

  Teuchos::RCP<Teuchos::Time> timer = Teuchos::rcp(new Teuchos::Time("1D Assembler global insert"));
  timer->start(true);

  // c a b
  for (LocalOrdinal i = 0; i < NumMyElements; i++) {
    if (MyGlobalElements[i] == indexBase) {
      // off-diagonal for first row
      Indices[0] = 1 + indexBase;

      // Diffusive term
      Values[0] = b * diff(xleft + MyGlobalElements[i] * dx) / (dx * dx);
      Values[0] = Values[0] - diff_prime(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);

      // Advective term
      Values[0] = Values[0] + adv(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);

      NumEntries = 1;

    } else if (MyGlobalElements[i] == NumGlobalElements + indexBase - 1) {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements - 2 + indexBase;

      // Diffusive term
      Values[0] = c * diff(xleft + MyGlobalElements[i] * dx) / (dx * dx);
      Values[0] = Values[0] + diff_prime(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);

      // Advective term
      Values[0] = Values[0] - adv(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);

      NumEntries = 1;

    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Indices[1] = MyGlobalElements[i] + 1;

      // Diffusive term
      Values[0] = c * diff(xleft + MyGlobalElements[i] * dx) / (dx * dx);
      Values[0] = Values[0] + diff_prime(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);
      Values[1] = b * diff(xleft + MyGlobalElements[i] * dx) / (dx * dx);
      Values[1] = Values[1] - diff_prime(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);

      // Advective term
      Values[0] = Values[0] - adv(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);
      Values[1] = Values[1] + adv(xleft + MyGlobalElements[i] * dx) / (2.0 * dx);

      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
    mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    // Diffusion
    Scalar diag_entry = a * diff(xleft + MyGlobalElements[i] * dx) / (dx * dx);
    // Reaction
    diag_entry = diag_entry + reac(xleft + MyGlobalElements[i] * dx);

    mtx->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(diag_entry));
  }

  timer->stop();

  timer = Teuchos::rcp(new Teuchos::Time("1D Assembler fillComplete"));
  timer->start(true);

  mtx->fillComplete();
  this->A_ = mtx;

  timer->stop();

  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

// =============================================  ADR2D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class ADR2DProblem : public Problem<Map, Matrix, MultiVector> {
 public:
  ADR2DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
    : Problem<Map, Matrix, MultiVector>(list, map) {}
  Teuchos::RCP<Matrix> BuildMatrix();

 private:
  // Function that defines the diffusion coefficient
  inline Scalar diff(Scalar x, Scalar y) { return 1.0 + 0.0 * x + 0.0 * y; };

  // Function that defines the first derivative in x-direction of the diffusion coefficient
  inline Scalar diff_primex(Scalar x, Scalar y) { return 0.0 + 0.0 * x + 0.0 * y; };

  // Function that defines the first derivative in x-direction of the diffusion coefficient
  inline Scalar diff_primey(Scalar x, Scalar y) { return 0.0 + 0.0 * x + 0.0 * y; };

  // Function that defines the advection coefficient in the x-direction
  inline Scalar advx(Scalar x, Scalar y) { return 10.0 + 0.0 * x + 0.0 * y; };

  // Function that defines the advection coefficient in the x-direction
  inline Scalar advy(Scalar x, Scalar y) { return 10.0 + 0.0 * x + 0.0 * y; };

  // Function taht defines the reaction coefficient
  inline Scalar reac(Scalar x, Scalar y) { return 0.0 + 0.0 * x + 0.0 * y; };
};

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> ADR2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal)-1);
  GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal)-1);
  double one       = 1.0;
  Scalar stretchx  = (Scalar)this->list_.get("stretchx", one);
  Scalar stretchy  = (Scalar)this->list_.get("stretchy", one);

  if (nx == -1 || ny == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)sqrt((double)n);
    ny              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny != n, std::logic_error, "You need to specify nx and ny.");
  }
  bool keepBCs = this->list_.get("keepBCs", false);

  // Diffusion stencil
  Scalar c1 = (Scalar)-one / (stretchx * stretchx);
  Scalar b1 = (Scalar)-one / (stretchx * stretchx);
  Scalar e1 = (Scalar)-one / (stretchy * stretchy);
  Scalar d1 = (Scalar)-one / (stretchy * stretchy);
  Scalar a1 = -(b1 + c1 + d1 + e1);

  // Advection stencil
  Scalar c2 = (Scalar)one / (2.0 * stretchx);
  Scalar b2 = (Scalar)one / (2.0 * stretchx);
  Scalar e2 = (Scalar)one / (2.0 * stretchy);
  Scalar d2 = (Scalar)one / (2.0 * stretchy);

  // this->A_ = Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, center, west, east, south, north, this->DirichletBC_, keepBCs);
  LocalOrdinal nnz = 5;

  Teuchos::RCP<Matrix> mtx = Galeri::Xpetra::MatrixTraits<Map, Matrix>::Build(this->Map_, nnz);

  LocalOrdinal numMyElements = (this->Map_)->getLocalNumElements();
  GlobalOrdinal indexBase    = (this->Map_)->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = (this->Map_)->getLocalElementList();

  GlobalOrdinal center, left, right, lower, upper;
  std::vector<Scalar> vals(nnz);
  std::vector<GlobalOrdinal> inds(nnz);

  //    e
  //  b a c
  //    d
  for (LocalOrdinal i = 0; i < numMyElements; ++i) {
    size_t n = 0;

    center = myGlobalElements[i] - indexBase;

    // Determine coordinates
    Scalar x1 = (Scalar)(center % nx) * stretchx;
    Scalar x2 = (Scalar)(std::floor(center / nx)) * stretchy;

    Galeri::Xpetra::GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);

    bool isDirichlet = (left == -1 && (this->DirichletBC_ & DIR_LEFT)) ||
                       (right == -1 && (this->DirichletBC_ & DIR_RIGHT)) ||
                       (lower == -1 && (this->DirichletBC_ & DIR_BOTTOM)) ||
                       (upper == -1 && (this->DirichletBC_ & DIR_TOP));

    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      inds[n]   = center;
      vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

    } else {
      // The Neumann b.c. are treated in a sane way. The Dirichlet b.c., however, are treated
      // insane when the option keepBCs=false. Speicifically, in this case we don't want to keep
      // Dirichlet b.c., but that would result in inconsistency between the map and the number of
      // degrees of freedom, plus the problem with GIDs. Therefore, we virtually expand domain by
      // one node in the direction of the Dirichlet b.c., and then assume that that node was
      // not kept. But we use an old GIDs. So yes, that's weird.

      if (left != -1) {
        inds[n]   = left;
        vals[n++] = b1 * diff(x1, x2) + b2 * diff_primex(x1, x2) - b2 * advx(x1, x2);
      }
      if (right != -1) {
        inds[n]   = right;
        vals[n++] = c1 * diff(x1, x2) - c2 * diff_primex(x1, x2) + c2 * advx(x1, x2);
      }
      if (lower != -1) {
        inds[n]   = lower;
        vals[n++] = d1 * diff(x1, x2) + d2 * diff_primey(x1, x2) - d2 * advy(x1, x2);
      }
      if (upper != -1) {
        inds[n]   = upper;
        vals[n++] = e1 * diff(x1, x2) - e2 * diff_primey(x1, x2) + e2 * advy(x1, x2);
      }

      // diagonal
      Scalar z = a1 * diff(x1, x2) + reac(x1, x2);
      if (Galeri::Xpetra::IsBoundary2d(center, nx, ny) && !isDirichlet) {
        // Neumann boundary unknown (diagonal = sum of all offdiagonal)
        z = Teuchos::ScalarTraits<Scalar>::zero();
        for (size_t j = 0; j < n; j++)
          z -= vals[j];
      }
      inds[n]   = center;
      vals[n++] = z;
    }

    for (size_t j = 0; j < n; j++)
      inds[j] += indexBase;

    Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
    Teuchos::ArrayView<Scalar> av(&vals[0], n);
    mtx->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  mtx->fillComplete();
  this->A_ = mtx;

  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

// =============================================  ADR3D  =============================================
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
class ADR3DProblem : public Problem<Map, Matrix, MultiVector> {
 public:
  ADR3DProblem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map)
    : Problem<Map, Matrix, MultiVector>(list, map) {}
  Teuchos::RCP<Matrix> BuildMatrix();

 private:
  // Function that defines the diffusion coefficient
  inline Scalar diff(Scalar x, Scalar y, Scalar z) { return 1.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function that defines the first derivative in x-direction of the diffusion coefficient
  inline Scalar diff_primex(Scalar x, Scalar y, Scalar z) { return 0.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function that defines the first derivative in x-direction of the diffusion coefficient
  inline Scalar diff_primey(Scalar x, Scalar y, Scalar z) { return 0.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function that defines the first derivative in x-direction of the diffusion coefficient
  inline Scalar diff_primez(Scalar x, Scalar y, Scalar z) { return 0.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function that defines the advection coefficient in the x-direction
  inline Scalar advx(Scalar x, Scalar y, Scalar z) { return 10.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function that defines the advection coefficient in the x-direction
  inline Scalar advy(Scalar x, Scalar y, Scalar z) { return 10.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function that defines the advection coefficient in the x-direction
  inline Scalar advz(Scalar x, Scalar y, Scalar z) { return 10.0 + 0.0 * x + 0.0 * y + 0.0 * z; };

  // Function taht defines the reaction coefficient
  inline Scalar reac(Scalar x, Scalar y, Scalar z) { return 0.0 + 0.0 * x + 0.0 * y + 0.0 * z; };
};

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
Teuchos::RCP<Matrix> ADR3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>::BuildMatrix() {
  GlobalOrdinal nx = this->list_.get("nx", (GlobalOrdinal)-1);
  GlobalOrdinal ny = this->list_.get("ny", (GlobalOrdinal)-1);
  GlobalOrdinal nz = this->list_.get("nz", (GlobalOrdinal)-1);
  double one       = 1.0;
  Scalar stretchx  = (Scalar)this->list_.get("stretchx", one);
  Scalar stretchy  = (Scalar)this->list_.get("stretchy", one);
  Scalar stretchz  = (Scalar)this->list_.get("stretchz", one);

  if (nx == -1 || ny == -1 || nz == -1) {
    GlobalOrdinal n = this->Map_->getGlobalNumElements();
    nx              = (GlobalOrdinal)Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny              = nx;
    nz              = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }
  bool keepBCs = this->list_.get("keepBCs", false);

  // Diffusion stencil
  Scalar c1 = (Scalar)-one / (stretchx * stretchx);
  Scalar b1 = (Scalar)-one / (stretchx * stretchx);
  Scalar g1 = (Scalar)-one / (stretchy * stretchy);
  Scalar f1 = (Scalar)-one / (stretchy * stretchy);
  Scalar e1 = (Scalar)-one / (stretchz * stretchz);
  Scalar d1 = (Scalar)-one / (stretchz * stretchz);
  Scalar a1 = -(c1 + b1 + g1 + f1 + e1 + d1);

  // Advection stencil
  Scalar c2 = (Scalar)one / (2.0 * stretchx);
  Scalar b2 = (Scalar)one / (2.0 * stretchx);
  Scalar g2 = (Scalar)one / (2.0 * stretchy);
  Scalar f2 = (Scalar)one / (2.0 * stretchy);
  Scalar e2 = (Scalar)one / (2.0 * stretchz);
  Scalar d2 = (Scalar)one / (2.0 * stretchz);

  // this->A_ = Cross3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(this->Map_, nx, ny, nz, center, left, right, front, back, down, up, this->DirichletBC_, keepBCs);
  LocalOrdinal nnz = 7;

  Teuchos::RCP<Matrix> mtx = Galeri::Xpetra::MatrixTraits<Map, Matrix>::Build(this->Map_, nnz);

  LocalOrdinal numMyElements = this->Map_->getLocalNumElements();
  GlobalOrdinal indexBase    = this->Map_->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = (this->Map_)->getLocalElementList();

  GlobalOrdinal center, left, right, bottom, top, front, back;
  std::vector<GlobalOrdinal> inds(nnz);
  std::vector<Scalar> vals(nnz);

  //    e
  //  b a c
  //    d
  // + f bottom and g top
  for (LocalOrdinal i = 0; i < numMyElements; ++i) {
    size_t n = 0;

    center = myGlobalElements[i] - indexBase;

    // Determine coordinates
    Scalar x3 = (Scalar)(std::floor(center / (nx * ny))) * stretchz;
    int plane = center % (nx * ny);
    Scalar x1 = (Scalar)(plane % nx) * stretchx;
    Scalar x2 = (Scalar)(std::floor(plane / nx)) * stretchy;

    Galeri::Xpetra::GetNeighboursCartesian3d(center, nx, ny, nz,
                                             left, right, front, back, bottom, top);

    bool isDirichlet = (left == -1 && (this->DirichletBC_ & DIR_LEFT)) ||
                       (right == -1 && (this->DirichletBC_ & DIR_RIGHT)) ||
                       (front == -1 && (this->DirichletBC_ & DIR_BOTTOM)) ||
                       (back == -1 && (this->DirichletBC_ & DIR_TOP)) ||
                       (front == -1 && (this->DirichletBC_ & DIR_FRONT)) ||
                       (back == -1 && (this->DirichletBC_ & DIR_BACK));

    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      inds[n]   = center;
      vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

    } else {
      // See comments about weird in Cross2D
      if (left != -1) {
        inds[n]   = left;
        vals[n++] = b1 * diff(x1, x2, x3) + b2 * diff_primex(x1, x2, x3) - b2 * advx(x1, x2, x3);
      }
      if (right != -1) {
        inds[n]   = right;
        vals[n++] = c1 * diff(x1, x2, x3) - c2 * diff_primex(x1, x2, x3) + c2 * advx(x1, x2, x3);
      }
      if (front != -1) {
        inds[n]   = front;
        vals[n++] = d1 * diff(x1, x2, x3) + d2 * diff_primey(x1, x2, x3) - d2 * advy(x1, x2, x3);
      }
      if (back != -1) {
        inds[n]   = back;
        vals[n++] = e1 * diff(x1, x2, x3) - e2 * diff_primey(x1, x2, x3) + e2 * advy(x1, x2, x3);
      }
      if (bottom != -1) {
        inds[n]   = bottom;
        vals[n++] = f1 * diff(x1, x2, x3) + f2 * diff_primez(x1, x2, x3) - f2 * advz(x1, x2, x3);
      }
      if (top != -1) {
        inds[n]   = top;
        vals[n++] = g1 * diff(x1, x2, x3) - g2 * diff_primez(x1, x2, x3) + g2 * advz(x1, x2, x3);
      }

      // diagonal
      Scalar z = a1 * diff(x1, x2, x3) + reac(x1, x2, x3);
      if (Galeri::Xpetra::IsBoundary3d(center, nx, ny, nz) && !isDirichlet) {
        // Neumann boundary unknown (diagonal = sum of all offdiagonal)
        z = Teuchos::ScalarTraits<Scalar>::zero();
        for (size_t j = 0; j < n; j++)
          z -= vals[j];
      }
      inds[n]   = center;
      vals[n++] = z;
    }

    for (size_t j = 0; j < n; j++)
      inds[j] += indexBase;

    Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
    Teuchos::ArrayView<Scalar> av(&vals[0], n);
    mtx->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  mtx->fillComplete();
  this->A_ = mtx;

  this->A_->setObjectLabel(this->getObjectLabel());
  return this->A_;
}

}  // namespace Xpetra

}  // namespace ADR

#endif  // CREATEADRMATRIX_HPP
