/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#include <fei_macros.hpp>

#include <test_utils/test_Factory_helper.hpp>

#include <test_utils/LibraryFactory.hpp>

#include <fei_Factory_Trilinos.hpp>

#include <snl_fei_Factory.hpp>

#include <fei_Vector_Impl.hpp>

#include <fei_Matrix_Impl.hpp>

#undef fei_file
#define fei_file "test_Factory_helper.cpp"
#include <fei_ErrMacros.hpp>

int test_Factory_helper::dyncastMatrix(fei::Matrix* matrix,
				       const char* libname)
{
  std::string sname(libname);

  if (sname == "TEST_LSC") {

    fei::Matrix_Impl<LinearSystemCore>* smatrix2 =
      dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matrix);
    if (smatrix2 == NULL) {
      fei::console_out() << "dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*> failed"<<FEI_ENDL;
      ERReturn(-1);
    }
  }

  if (sname == "Aztec") {
#ifdef HAVE_FEI_AZTECOO
    fei::Matrix_Impl<LinearSystemCore>* smatrix =
      dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matrix);
    if (smatrix == NULL) {
      fei::console_out() << "dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*> failed"<<FEI_ENDL;
      ERReturn(-1);
    }
#else
    fei::console_out() << "libname==Aztec but HAVE_FEI_AZTECOO not defined."<<FEI_ENDL;
    ERReturn(-1);
#endif
  }

  if (sname == "Trilinos") {
#ifdef HAVE_FEI_EPETRA
    fei::Matrix_Impl<Epetra_CrsMatrix>* smatrix =
      dynamic_cast<fei::Matrix_Impl<Epetra_CrsMatrix>*>(matrix);
    if (smatrix == NULL) {
      fei::console_out() << "dynamic_cast<fei::Matrix_Impl<Epetra_CrsMatrix>*> failed"<<FEI_ENDL;
      ERReturn(-1);
    }
#else
    fei::console_out() << "libname==Trilinos but HAVE_FEI_EPETRA not defined."<<FEI_ENDL;
    ERReturn(-1);
#endif
  }

  return(0);
}

int test_Factory_helper::dyncastVector(fei::Vector* vector,
				       const char* libname)
{
  std::string sname(libname);
  if (sname == "TEST_LSC") {
    fei::Vector_Impl<LinearSystemCore>* svector =
      dynamic_cast<fei::Vector_Impl<LinearSystemCore>*>(vector);
    if (svector == NULL) {
      fei::console_out() << "dynamic_cast<fei::Vector_Impl<LinearSystemCore>*> failed"<<FEI_ENDL;
      ERReturn(-1);
    }
  }

  if (sname == "Aztec") {
#ifdef HAVE_FEI_AZTECOO
    fei::Vector_Impl<LinearSystemCore>* svector =
      dynamic_cast<fei::Vector_Impl<LinearSystemCore>*>(vector);
    if (svector == NULL) {
      fei::console_out() << "dynamic_cast<fei::Vector_Impl<LinearSystemCore>*> failed"<<FEI_ENDL;
      ERReturn(-1);
    }
#else
    fei::console_out() << "libname==Aztec but HAVE_FEI_AZTECOO not defined."<<FEI_ENDL;
    ERReturn(-1);
#endif
  }

  if (sname == "Trilinos") {
#ifdef HAVE_FEI_EPETRA
    fei::Vector_Impl<Epetra_MultiVector>* svector =
      dynamic_cast<fei::Vector_Impl<Epetra_MultiVector>*>(vector);
    if (svector == NULL) {
      fei::console_out() << "dynamic_cast<fei::Vector_Impl<Epetra_MultiVector>*> failed"<<FEI_ENDL;
      ERReturn(-1);
    }
#else
    fei::console_out() << "libname==Trilinos but HAVE_FEI_EPETRA not defined."<<FEI_ENDL;
    ERReturn(-1);
#endif
  }

  return(0);
}
