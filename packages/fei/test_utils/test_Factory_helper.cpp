/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
