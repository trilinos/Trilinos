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

#ifndef __TEKO_CRS_MATRIX_GALLERY_TPETRA_HPP
#define __TEKO_CRS_MATRIX_GALLERY_TPETRA_HPP

#include "Teuchos_RCPDecl.hpp"
#include "Trilinos_Util_CommandLineParser.h"
#include <string>
#include <vector>

#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_Vector_decl.hpp"

#include "Teko_ConfigDefs.hpp"

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

class CrsMatrixGallery {
public:
  CrsMatrixGallery(const std::string &name);
  ~CrsMatrixGallery();

  //! Sets a gallery options using an interger value.
  int Set(const std::string &parameter, const int value);

  //! Returns a pointer to the CrsMatrix.
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> GetMatrix();

  void CreateMap();

  //! Creates the CrdMatrix.
  void CreateMatrix();

  // Creates a diagonal matrix. Elements on the diagonal are called `a'.
  void CreateMatrixDiag();
  void CreateMatrixLaplace2d();
  void CreateMatrixRecirc2d();

  // returns the neighbors of a given node. The node is supposed to be on
  // a 2D Cartesian grid
  void GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                                int &left, int &right, int &lower, int &upper);

  // put to NULL or default values all internal data
  void ZeroOutData();

  void SetupCartesianGrid2D();

  void CreateMatrixCrossStencil2d();
  void CreateMatrixCrossStencil2dVector();

  //@}

  // ======================== //
  // I N T E R N A L  D A T A //
  // ======================== //

  Teuchos::RCP<const Teuchos::Comm<int>> comm_;

  // matrix and vectors (scalar)
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> matrix_;
  Teuchos::RCP<Tpetra::Map<LO, GO>> map_;

  // information about the problem to generate
  std::string name_;
  GO NumGlobalElements_;
  LO NumMyElements_;

  // parameters
  int nx_, ny_;

  double lx_, ly_;

  Teuchos::RCP<Tpetra::Vector<ST, LO, GO, NT>> VectorA_, VectorB_, VectorC_,
      VectorD_, VectorE_;

  double a_, b_, c_, d_, e_;
  double conv_, diff_;

  // others
  std::string ErrorMsg;
  std::string OutputMsg;
  bool verbose_;
};

#endif // __TEKO_CRS_MATRIX_GALLERY_TPETRA_HPP
