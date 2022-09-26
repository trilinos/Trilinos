// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_TestForException.hpp"
#include "Trilinos_Util.h"
#include <string>

#include "CrsMatrixGalleryTpetra.hpp"
#include "Trilinos_Util_CommandLineParser.h"

#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_Map_def.hpp"
#include "Tpetra_MultiVector_def.hpp"
#include "Tpetra_Vector_def.hpp"

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

const double UNDEF = -99999.87;
const bool Scaling = false;

// ================================================ ====== ==== ==== == =
CrsMatrixGallery::CrsMatrixGallery(const string &name) : name_(name) {
  comm_ = Tpetra::getDefaultComm();
  ZeroOutData();
  // verbosity level (now always false)
  // if( comm_->MyPID()==0 ) verbose_ = true;
  verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [CrsMatrixGalleryTpetra]: ";
  OutputMsg = "CrsMatrixGalleryTpetra: ";
}

// ================================================ ====== ==== ==== == =
CrsMatrixGallery::~CrsMatrixGallery() {
  // put to default values
  ZeroOutData();
}

// ================================================ ====== ==== ==== == =
int CrsMatrixGallery::Set(const std::string &parameter, const int value) {
  if (parameter == "nx") {

    if (value <= 0) {
      cerr << ErrorMsg << "nx must be greater than 0\n";
      return -1;
    }

    nx_ = value;
    return 0;

  } else if (parameter == "ny") {

    if (value <= 0) {
      cerr << ErrorMsg << "ny must be greater than 0\n";
      return -1;
    }

    ny_ = value;
    return 0;
  }

  cerr << ErrorMsg << "input string (" << parameter << ") not valid\n";
  return -2;
}

// ================================================ ====== ==== ==== == =
RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> CrsMatrixGallery::GetMatrix() {
  if (matrix_ == Teuchos::null)
    CreateMatrix();
  return matrix_;
}

void CrsMatrixGallery::CreateMap(void) {
  // first get the problem size. For some problems. the user can
  // specify the problem size using different parameters (e.g.,
  // nx and ny for a 2D Laplace problem). I need the internal
  // variable NumGlobalElements_ properly set before continuing.

  if (name_ == "diag") {
    if (NumGlobalElements_ <= 0) {
      if (nx_ > 0)
        NumGlobalElements_ = nx_;
      else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, ErrorMsg
                                              << "problem size not correct ("
                                              << NumGlobalElements_ << ")");
      }
    }
  } else if (name_ == "laplace_2d" || name_ == "recirc_2d") {

    if (NumGlobalElements_ <= 0) {
      if (nx_ > 0 && ny_ > 0)
        NumGlobalElements_ = nx_ * ny_;
      else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true, ErrorMsg << "Problem size not correct (" << NumGlobalElements_
                           << ")\n"
                           << ErrorMsg << "It should be a perfect square");
      }
    }

    if (verbose_) {
      cout << OutputMsg << "nx = " << nx_ << ", ny = " << ny_ << endl;
    }

  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        ErrorMsg << "matrix name is incorrect or not set (" << name_ << ")");
  }

  map_ = rcp(new Tpetra::Map<LO, GO>(NumGlobalElements_, GO{0}, comm_));

  // local number of rows
  NumMyElements_ = map_->getLocalNumElements();
}

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::CreateMatrix(void) {
  if (map_ == Teuchos::null)
    CreateMap();

  if (name_ == "diag")
    CreateMatrixDiag();
  else if (name_ == "laplace_2d")
    CreateMatrixLaplace2d();
  else if (name_ == "recirc_2d")
    CreateMatrixRecirc2d();
  else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        ErrorMsg << "matrix name is incorrect or not set (" << name_ << ")");
  }
}

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::CreateMatrixDiag(void) {

  // default value if not otherwise specified
  if (a_ == UNDEF)
    a_ = 1;

  if (verbose_) {
    cout << OutputMsg << "Creating matrix `diag'...\n";
    cout << OutputMsg << "Diagonal element = " << a_ << endl;
  }

  matrix_ = rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map_, 1));

  for (int i = 0; i < NumMyElements_; ++i) {
    auto index = map_->getGlobalElement(i);

    matrix_->insertGlobalValues(index, Teuchos::tuple<GO>(1),
                                Teuchos::tuple<ST>(a_));
  }

  matrix_->fillComplete();
}

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::CreateMatrixCrossStencil2d(void) {
  // default values if not otherwise specified

  if (a_ == UNDEF)
    a_ = 4;
  if (b_ == UNDEF)
    b_ = 1;
  if (c_ == UNDEF)
    c_ = 1;
  if (d_ == UNDEF)
    d_ = 1;
  if (e_ == UNDEF)
    e_ = 1;

  if (verbose_) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
    cout << OutputMsg << "with values: a=" << a_ << ", b=" << b_ << ", c=" << c_
         << ", d=" << d_ << ", e=" << e_ << endl;
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map_, 5));

  // Add  rows one-at-a-time

  double Values[4], diag;
  GO Indices[4];

  //    e
  //  b a c
  //    d
  const auto myGlobalElems = map_->getLocalElementList();

  for (int i = 0; i < NumMyElements_; ++i) {
    int NumEntries = 0;
    GetNeighboursCartesian2d(myGlobalElems[i], nx_, ny_, left, right, lower,
                             upper);
    if (left != -1) {
      Indices[NumEntries] = left;
      Values[NumEntries] = b_;
      ++NumEntries;
    }
    if (right != -1) {
      Indices[NumEntries] = right;
      Values[NumEntries] = c_;
      ++NumEntries;
    }
    if (lower != -1) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d_;
      ++NumEntries;
    }
    if (upper != -1) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e_;
      ++NumEntries;
    }
    // put the off-diagonal entries
    matrix_->insertGlobalValues(myGlobalElems[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    diag = a_;

    matrix_->insertGlobalValues(myGlobalElems[i], 1, &diag, &myGlobalElems[i]);
  }
  matrix_->fillComplete();
}

void CrsMatrixGallery::CreateMatrixLaplace2d() {

  SetupCartesianGrid2D();

  double hx = lx_ / (nx_ + 1);
  double hy = ly_ / (ny_ + 1);

  if (verbose_) {
    cout << OutputMsg << "Creating matrix `laplace_2d'...\n";
  }

  if (Scaling) {
    if (Scaling) {
      cout << OutputMsg << "hx = " << hx << ", hy = " << hy << endl;
    }
    a_ = 2.0 / (hx * hx) + 2.0 / (hy * hy);
    b_ = -1.0 / (hx * hx);
    c_ = -1.0 / (hx * hx);
    d_ = -1.0 / (hy * hy);
    e_ = -1.0 / (hy * hy);
  } else {
    a_ = 4;
    b_ = -1;
    c_ = -1;
    d_ = -1;
    e_ = -1;
  }

  CreateMatrixCrossStencil2d();
}

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::CreateMatrixRecirc2d(void) {

  // default values if not specified otherwise

  if (conv_ == UNDEF)
    conv_ = 1;
  if (diff_ == UNDEF)
    diff_ = 1e-5;

  if (verbose_) {
    cout << OutputMsg << "Creating matrix `recirc_2d'...\n";
    cout << OutputMsg << "with convection = " << conv_
         << " and diffusion = " << diff_ << endl;
  }

  SetupCartesianGrid2D();

  if (VectorA_ == Teuchos::null)
    VectorA_ = rcp(new Tpetra::Vector<ST, LO, GO, NT>(map_));
  if (VectorB_ == Teuchos::null)
    VectorB_ = rcp(new Tpetra::Vector<ST, LO, GO, NT>(map_));
  if (VectorC_ == Teuchos::null)
    VectorC_ = rcp(new Tpetra::Vector<ST, LO, GO, NT>(map_));
  if (VectorD_ == Teuchos::null)
    VectorD_ = rcp(new Tpetra::Vector<ST, LO, GO, NT>(map_));
  if (VectorE_ == Teuchos::null)
    VectorE_ = rcp(new Tpetra::Vector<ST, LO, GO, NT>(map_));

  VectorA_->putScalar(0.0);
  VectorB_->putScalar(0.0);
  VectorC_->putScalar(0.0);
  VectorD_->putScalar(0.0);
  VectorE_->putScalar(0.0);

  double hx = lx_ / (nx_ + 1);
  double hy = ly_ / (ny_ + 1);

  // int_type *&MyGlobalElements = MyGlobalElementsPtr<int_type>();
  const auto myGlobalElems = map_->getLocalElementList();

  {
    auto vecA_2d = VectorA_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto vecA_1d = Kokkos::subview(vecA_2d, Kokkos::ALL, 0);

    auto vecB_2d = VectorB_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto vecB_1d = Kokkos::subview(vecB_2d, Kokkos::ALL, 0);

    auto vecC_2d = VectorC_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto vecC_1d = Kokkos::subview(vecC_2d, Kokkos::ALL, 0);

    auto vecD_2d = VectorD_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto vecD_1d = Kokkos::subview(vecD_2d, Kokkos::ALL, 0);

    auto vecE_2d = VectorE_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto vecE_1d = Kokkos::subview(vecE_2d, Kokkos::ALL, 0);

    for (int i = 0; i < NumMyElements_; ++i) {
      int ix, iy;
      ix = (myGlobalElems[i]) % nx_;
      iy = (myGlobalElems[i] - ix) / nx_;
      double x = hx * (ix + 1);
      double y = hy * (iy + 1);
      double ConvX = conv_ * 4 * x * (x - 1.) * (1. - 2 * y) / hx;
      double ConvY = -conv_ * 4 * y * (y - 1.) * (1. - 2 * x) / hy;

      // convection part

      if (ConvX < 0) {
        vecC_1d(i) += ConvX;
        vecA_1d(i) -= ConvX;
      } else {
        vecB_1d(i) -= ConvX;
        vecA_1d(i) += ConvX;
      }

      if (ConvY < 0) {
        vecE_1d(i) += ConvY;
        vecA_1d(i) -= ConvY;
      } else {
        vecD_1d(i) -= ConvY;
        vecA_1d(i) += ConvY;
      }

      // add diffusion part
      vecA_1d(i) += diff_ * 2. / (hx * hx) + diff_ * 2. / (hy * hy);
      vecB_1d(i) -= diff_ / (hx * hx);
      vecC_1d(i) -= diff_ / (hx * hx);
      vecD_1d(i) -= diff_ / (hy * hy);
      vecE_1d(i) -= diff_ / (hy * hy);
    }
  }

  CreateMatrixCrossStencil2dVector();
}

void CrsMatrixGallery::CreateMatrixCrossStencil2dVector() {

  if (verbose_) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
  }

  SetupCartesianGrid2D();

  int left, right, lower, upper;

  matrix_ = rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map_, 5));

  // Add  rows one-at-a-time

  double Values[4], diag;
  GO Indices[4];

  //    e
  //  b a c
  //    d
  // int_type*& MyGlobalElements = MyGlobalElementsPtr<int_type>();
  const auto myGlobalElems = map_->getLocalElementList();

  auto vecA_2d = VectorA_->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto vecA_1d = Kokkos::subview(vecA_2d, Kokkos::ALL, 0);

  auto vecB_2d = VectorB_->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto vecB_1d = Kokkos::subview(vecB_2d, Kokkos::ALL, 0);

  auto vecC_2d = VectorC_->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto vecC_1d = Kokkos::subview(vecC_2d, Kokkos::ALL, 0);

  auto vecD_2d = VectorD_->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto vecD_1d = Kokkos::subview(vecD_2d, Kokkos::ALL, 0);

  auto vecE_2d = VectorE_->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto vecE_1d = Kokkos::subview(vecE_2d, Kokkos::ALL, 0);

  for (int i = 0; i < NumMyElements_; ++i) {
    int NumEntries = 0;
    GetNeighboursCartesian2d(myGlobalElems[i], nx_, ny_, left, right, lower,
                             upper);
    if (left != -1) {
      Indices[NumEntries] = left;
      Values[NumEntries] = vecB_1d(i);
      ++NumEntries;
    }
    if (right != -1) {
      Indices[NumEntries] = right;
      Values[NumEntries] = vecC_1d(i);
      ++NumEntries;
    }
    if (lower != -1) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = vecD_1d(i);
      ++NumEntries;
    }
    if (upper != -1) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = vecE_1d(i);
      ++NumEntries;
    }
    // put the off-diagonal entries
    matrix_->insertGlobalValues(myGlobalElems[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    diag = vecA_1d(i);

    matrix_->insertGlobalValues(myGlobalElems[i], 1, &diag, &myGlobalElems[i]);
  }

  matrix_->fillComplete();
}

// ================================================ ====== ==== ==== == =

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::GetNeighboursCartesian2d(const int i, const int nx,
                                                const int ny, int &left,
                                                int &right, int &lower,
                                                int &upper) {

  int ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)
    left = -1;
  else
    left = i - 1;
  if (ix == nx - 1)
    right = -1;
  else
    right = i + 1;
  if (iy == 0)
    lower = -1;
  else
    lower = i - nx;
  if (iy == ny - 1)
    upper = -1;
  else
    upper = i + nx;

  return;
}

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::ZeroOutData() {
  NumGlobalElements_ = -1;
  nx_ = -1;
  ny_ = -1;

  lx_ = 1.0;
  ly_ = 1.0;

  a_ = UNDEF, b_ = UNDEF, c_ = UNDEF, d_ = UNDEF, e_ = UNDEF;

  conv_ = UNDEF;
  diff_ = UNDEF;

  VectorA_ = Teuchos::null;
  VectorB_ = Teuchos::null;
  VectorC_ = Teuchos::null;
  VectorD_ = Teuchos::null;
  VectorE_ = Teuchos::null;

  map_ = Teuchos::null;
  matrix_ = Teuchos::null;
}

ostream &operator<<(ostream &os, const CrsMatrixGallery &G) {

  bool verbose = (G.comm_->getRank() == 0);

  if (verbose) {

    os << " * Solving problem " << G.name_ << endl;
    os << " * Number of global elements : " << G.NumGlobalElements_ << endl;

    // CRS stuff
    if (G.matrix_ != Teuchos::null) {
      os << " * the matrix has been created " << endl;
      os << " * Matrix->getDomainMap()->getGlobalNumElements() = "
         << G.matrix_->getDomainMap()->getGlobalNumElements() << endl;
    }
  }

  return os;
}

// ================================================ ====== ==== ==== == =
void CrsMatrixGallery::SetupCartesianGrid2D() {
  // needs a square number of nodes or
  // nx and ny set
  if (nx_ == -1 || ny_ == -1) {
    nx_ = (int)sqrt((double)NumGlobalElements_);
    ny_ = nx_;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        nx_ * ny_ != NumGlobalElements_,
        ErrorMsg << "The number of global elements must be a perfect square\n"
                 << ErrorMsg << "otherwise set nx and ny.\n"
                 << ErrorMsg
                 << "(now NumGlobalElements = " << NumGlobalElements_ << ")\n");
  }
}