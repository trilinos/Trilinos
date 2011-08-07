// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// System includes
#include <algorithm>

// Epetra includes
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
// #include "Epetra_FastCrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
//#include "Epetra_MsrMatrix.h"
#include "Epetra_VbrRowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_JadMatrix.h"

// Local includes
#include "PyTrilinos_config.h"
#include "PyTrilinos_Util.h"
#include "Epetra_PyUtil.h"
#include "PythonException.h"
#include "swigpyrun.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS
PyObject *
PyTrilinos::convertEpetraMultiVectorToPython(const Teuchos::RCP< Epetra_MultiVector > *emv)
{
  // SWIG initialization
  static swig_type_info * swig_ENMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< PyTrilinos::Epetra_NumPyMultiVector >*");
  //
  // Convert to PyTrilinos::Epetra_NumPyMultiVector
  const Teuchos::RCP< PyTrilinos::Epetra_NumPyMultiVector > *enmv = new
    Teuchos::RCP< PyTrilinos::Epetra_NumPyMultiVector >
    (new PyTrilinos::Epetra_NumPyMultiVector(View, **emv));
  return SWIG_NewPointerObj((void*)enmv, swig_ENMV_ptr, 1);
}
PyObject *
PyTrilinos::convertEpetraMultiVectorToPython(const Teuchos::RCP< const Epetra_MultiVector > *cemv)
{
  // SWIG initialization
  static swig_type_info * swig_ENMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< PyTrilinos::Epetra_NumPyMultiVector >*");
  //
  // Convert to PyTrilinos::Epetra_NumPyMultiVector
  const Teuchos::RCP< PyTrilinos::Epetra_NumPyMultiVector > *enmv = new
    Teuchos::RCP< PyTrilinos::Epetra_NumPyMultiVector >
    (new PyTrilinos::Epetra_NumPyMultiVector(View, **cemv));
  return SWIG_NewPointerObj((void*)enmv, swig_ENMV_ptr, 1);
}

#else

PyObject *
PyTrilinos::convertEpetraMultiVectorToPython(const Epetra_MultiVector *emv)
{
  // SWIG initialization
  static swig_type_info * swig_ENMV_ptr = SWIG_TypeQuery("PyTrilinos::Epetra_NumPyMultiVector *");
  //
  // Convert to PyTrilinos::Epetra_NumPyMultiVector
  const PyTrilinos::Epetra_NumPyMultiVector *enmv =
    new PyTrilinos::Epetra_NumPyMultiVector(View, *emv);
  return SWIG_NewPointerObj((void*)enmv, swig_ENMV_ptr, 1);
}
#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS
PyObject *
PyTrilinos::convertEpetraVectorToPython(const Teuchos::RCP< Epetra_Vector > *ev)
{
  // SWIG initialization
  static swig_type_info * swig_ENV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< PyTrilinos::Epetra_NumPyVector >*");
  //
  // Convert to PyTrilinos::Epetra_NumPyVector
  const Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *env = new
    Teuchos::RCP< PyTrilinos::Epetra_NumPyVector >(new PyTrilinos::Epetra_NumPyVector(View, **ev));
  return SWIG_NewPointerObj((void*)env, swig_ENV_ptr, 1);
}

PyObject *
PyTrilinos::convertEpetraVectorToPython(const Teuchos::RCP< const Epetra_Vector > *cev)
{
  // SWIG initialization
  static swig_type_info * swig_ENV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< PyTrilinos::Epetra_NumPyVector >*");
  //
  // Convert to PyTrilinos::Epetra_NumPyVector
  const Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *env = new
    Teuchos::RCP< PyTrilinos::Epetra_NumPyVector >(new PyTrilinos::Epetra_NumPyVector(View, **cev));
  return SWIG_NewPointerObj((void*)env, swig_ENV_ptr, 1);
}

#else

PyObject *
PyTrilinos::convertEpetraVectorToPython(const Epetra_Vector *emv)
{
  // SWIG initialization
  static swig_type_info * swig_ENV_ptr = SWIG_TypeQuery("PyTrilinos::Epetra_NumPyVector *");
  //
  const PyTrilinos::Epetra_NumPyVector *enmv = new PyTrilinos::Epetra_NumPyVector(View, *emv);
  return SWIG_NewPointerObj((void*)enmv, swig_ENV_ptr, 1);
}

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS

PyObject *
PyTrilinos::convertEpetraOperatorToPython(const Teuchos::RCP< Epetra_Operator > *eo)
{
  // SWIG initialization
  static swig_type_info *swig_EO_ptr   = SWIG_TypeQuery("Teuchos::RCP< Epetra_Operator        >*");
  //static swig_type_info *swig_EFCO_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_FastCrsOperator >*");
  static swig_type_info *swig_EIO_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_InvOperator     >*");
  static swig_type_info *swig_ERM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_RowMatrix       >*");
  static swig_type_info *swig_EBRM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_BasicRowMatrix  >*");
  static swig_type_info *swig_ECM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_CrsMatrix       >*");
  //static swig_type_info *swig_EMM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_MsrMatrix       >*");
  static swig_type_info *swig_EVM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrMatrix       >*");
  static swig_type_info *swig_EVRM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrRowMatrix    >*");
  static swig_type_info *swig_EFVM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_FEVbrMatrix     >*");
  static swig_type_info *swig_EFCM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_FECrsMatrix     >*");
  static swig_type_info *swig_EJM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_JadMatrix       >*");
  //
  // Attempt to downcast to Epetra_VbrRowMatrix
  Teuchos::RCP< Epetra_VbrRowMatrix > *evrm = new
    Teuchos::RCP< Epetra_VbrRowMatrix >(Teuchos::rcp_dynamic_cast< Epetra_VbrRowMatrix >(*eo));
  if (evrm->is_null()) delete evrm;
  else return SWIG_NewPointerObj((void*)evrm, swig_EVRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FEVbrMatrix
  Teuchos::RCP< Epetra_FEVbrMatrix > *efvm = new
    Teuchos::RCP< Epetra_FEVbrMatrix >(Teuchos::rcp_dynamic_cast< Epetra_FEVbrMatrix >(*eo));
  if (efvm->is_null()) delete efvm;
  else return SWIG_NewPointerObj((void*)efvm, swig_EFVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FECrsMatrix
  Teuchos::RCP< Epetra_FECrsMatrix > *efcm = new
    Teuchos::RCP< Epetra_FECrsMatrix >(Teuchos::rcp_dynamic_cast< Epetra_FECrsMatrix >(*eo));
  if (efcm->is_null()) delete efcm;
  else return SWIG_NewPointerObj((void*)efcm, swig_EFCM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_JadMatrix
  Teuchos::RCP< Epetra_JadMatrix > *ejm = new
    Teuchos::RCP< Epetra_JadMatrix >(Teuchos::rcp_dynamic_cast< Epetra_JadMatrix >(*eo));
  if (ejm->is_null()) delete ejm;
  else return SWIG_NewPointerObj((void*)ejm, swig_EJM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_BasicRowMatrix
  Teuchos::RCP< Epetra_BasicRowMatrix > *ebrm = new
    Teuchos::RCP< Epetra_BasicRowMatrix >(Teuchos::rcp_dynamic_cast< Epetra_BasicRowMatrix >(*eo));
  if (ebrm->is_null()) delete ebrm;
  else return SWIG_NewPointerObj((void*)ebrm, swig_EBRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_CrsMatrix
  Teuchos::RCP< Epetra_CrsMatrix > *ecm = new
    Teuchos::RCP< Epetra_CrsMatrix >(Teuchos::rcp_dynamic_cast< Epetra_CrsMatrix >(*eo));
  if (ecm->is_null()) delete ecm;
  else return SWIG_NewPointerObj((void*)ecm, swig_ECM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_MsrMatrix
  // Teuchos::RCP< Epetra_MsrMatrix > *emm = new
  //   Teuchos::RCP< Epetra_MsrMatrix >(Teuchos::rcp_dynamic_cast< Epetra_MsrMatrix >(*eo));
  // if (emm->is_null()) delete emm;
  // else return SWIG_NewPointerObj((void*)emm, swig_EMM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_VbrMatrix
  Teuchos::RCP< Epetra_VbrMatrix > *evm = new
    Teuchos::RCP< Epetra_VbrMatrix >(Teuchos::rcp_dynamic_cast< Epetra_VbrMatrix >(*eo));
  if (evm->is_null()) delete evm;
  else return SWIG_NewPointerObj((void*)evm, swig_EVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_RowMatrix
  Teuchos::RCP< Epetra_RowMatrix > *erm = new
    Teuchos::RCP< Epetra_RowMatrix >(Teuchos::rcp_dynamic_cast< Epetra_RowMatrix >(*eo));
  if (erm->is_null()) delete erm;
  else return SWIG_NewPointerObj((void*)erm, swig_ERM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_InvOperator
  Teuchos::RCP< Epetra_InvOperator > *eio = new
    Teuchos::RCP< Epetra_InvOperator >(Teuchos::rcp_dynamic_cast< Epetra_InvOperator >(*eo));
  if (eio->is_null()) delete eio;
  else return SWIG_NewPointerObj((void*)eio, swig_EIO_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FastCrsOperator
  // Teuchos::RCP< Epetra_FastCrsOperator > *efco = new
  //   Teuchos::RCP< Epetra_FastCrsOperator >(Teuchos::rcp_dynamic_cast< Epetra_FastCrsOperator >(*eo));
  // if (efco->is_null()) delete efco;
  // else return SWIG_NewPointerObj((void*)efco, swig_EFCO_ptr, SWIG_POINTER_OWN);
  //
  // If downcast is not possible, return python Epetra_Operator
  // **NOTE** Always make a copy!
  Teuchos::RCP< Epetra_Operator > *eocopy = new Teuchos::RCP< Epetra_Operator >(*eo);
  return SWIG_NewPointerObj((void*)eocopy, swig_EO_ptr, SWIG_POINTER_OWN);
}

PyObject *
PyTrilinos::convertEpetraOperatorToPython(const Teuchos::RCP< const Epetra_Operator > *ceo)
{
  // SWIG initialization
  static swig_type_info *swig_EO_ptr   = SWIG_TypeQuery("Teuchos::RCP< Epetra_Operator        >*");
  //static swig_type_info *swig_EFCO_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_FastCrsOperator >*");
  static swig_type_info *swig_EIO_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_InvOperator     >*");
  static swig_type_info *swig_ERM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_RowMatrix       >*");
  static swig_type_info *swig_EBRM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_BasicRowMatrix  >*");
  static swig_type_info *swig_ECM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_CrsMatrix       >*");
  //static swig_type_info *swig_EMM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_MsrMatrix       >*");
  static swig_type_info *swig_EVM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrMatrix       >*");
  static swig_type_info *swig_EVRM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrRowMatrix    >*");
  static swig_type_info *swig_EFVM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_FEVbrMatrix     >*");
  static swig_type_info *swig_EFCM_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_FECrsMatrix     >*");
  static swig_type_info *swig_EJM_ptr  = SWIG_TypeQuery("Teuchos::RCP< Epetra_JadMatrix       >*");
  //
  // Cast const-ness away
  Teuchos::RCP< Epetra_Operator > *eo = new
    Teuchos::RCP< Epetra_Operator >(Teuchos::rcp_const_cast< Epetra_Operator >(*ceo));
  //
  // Attempt to downcast to Epetra_VbrMatrix
  Teuchos::RCP< const Epetra_VbrRowMatrix > *evrm = new
    Teuchos::RCP< const Epetra_VbrRowMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_VbrRowMatrix >(*eo));
  if (evrm->is_null()) delete evrm;
  else return SWIG_NewPointerObj((void*)evrm, swig_EVRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FEVbrMatrix
  Teuchos::RCP< const Epetra_FEVbrMatrix > *efvm = new
    Teuchos::RCP< const Epetra_FEVbrMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_FEVbrMatrix >(*eo));
  if (efvm->is_null()) delete efvm;
  else return SWIG_NewPointerObj((void*)efvm, swig_EFVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FECrsMatrix
  Teuchos::RCP< const Epetra_FECrsMatrix > *efcm = new
    Teuchos::RCP< const Epetra_FECrsMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_FECrsMatrix >(*eo));
  if (efcm->is_null()) delete efcm;
  else return SWIG_NewPointerObj((void*)efcm, swig_EFCM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_JadMatrix
  Teuchos::RCP< const Epetra_JadMatrix > *ejm = new
    Teuchos::RCP< const Epetra_JadMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_JadMatrix >(*eo));
  if (ejm->is_null()) delete ejm;
  else return SWIG_NewPointerObj((void*)ejm, swig_EJM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_BasicRowMatrix
  Teuchos::RCP< const Epetra_BasicRowMatrix > *ebrm = new
    Teuchos::RCP< const Epetra_BasicRowMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_BasicRowMatrix >(*eo));
  if (ebrm->is_null()) delete ebrm;
  else return SWIG_NewPointerObj((void*)ebrm, swig_EBRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_CrsMatrix
  Teuchos::RCP< const Epetra_CrsMatrix > *ecm = new
    Teuchos::RCP< const Epetra_CrsMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_CrsMatrix >(*eo));
  if (ecm->is_null()) delete ecm;
  else return SWIG_NewPointerObj((void*)ecm, swig_ECM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_MsrMatrix
  // Teuchos::RCP< const Epetra_MsrMatrix > *emm = new
  //   Teuchos::RCP< const Epetra_MsrMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_MsrMatrix >(*eo));
  // if (emm->is_null()) delete emm;
  // else return SWIG_NewPointerObj((void*)emm, swig_EMM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_VbrMatrix
  Teuchos::RCP< const Epetra_VbrMatrix > *evm = new
    Teuchos::RCP< const Epetra_VbrMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_VbrMatrix >(*eo));
  if (evm->is_null()) delete evm;
  else return SWIG_NewPointerObj((void*)evm, swig_EVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_RowMatrix
  Teuchos::RCP< const Epetra_RowMatrix > *erm = new
    Teuchos::RCP< const Epetra_RowMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_RowMatrix >(*eo));
  if (erm->is_null()) delete erm;
  else return SWIG_NewPointerObj((void*)erm, swig_ERM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_InvOperator
  Teuchos::RCP< const Epetra_InvOperator > *eio = new
    Teuchos::RCP< const Epetra_InvOperator >(Teuchos::rcp_dynamic_cast< const Epetra_InvOperator >(*eo));
  if (eio->is_null()) delete eio;
  else return SWIG_NewPointerObj((void*)eio, swig_EIO_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FastCrsOperator
  // Teuchos::RCP< const Epetra_FastCrsOperator > *efco = new
  //   Teuchos::RCP< const Epetra_FastCrsOperator >(Teuchos::rcp_dynamic_cast< const Epetra_FastCrsOperator >(*eo));
  // if (efco->is_null()) delete efco;
  // else return SWIG_NewPointerObj((void*)efco, swig_EFCO_ptr, SWIG_POINTER_OWN);
  //
  // If downcast is not possible, return python Epetra_Operator
  return SWIG_NewPointerObj((void*)eo, swig_EO_ptr, SWIG_POINTER_OWN);
}

#else

PyObject *
PyTrilinos::convertEpetraOperatorToPython(const Epetra_Operator * eo, int cnvt_flags)
{
  // SWIG initialization
  static swig_type_info * swig_EO_ptr   = SWIG_TypeQuery("Epetra_Operator        *");
  //static swig_type_info * swig_EFCO_ptr = SWIG_TypeQuery("Epetra_FastCrsOperator *");
  static swig_type_info * swig_EIO_ptr  = SWIG_TypeQuery("Epetra_InvOperator     *");
  static swig_type_info * swig_ERM_ptr  = SWIG_TypeQuery("Epetra_RowMatrix       *");
  static swig_type_info * swig_EBRM_ptr = SWIG_TypeQuery("Epetra_BasicRowMatrix  *");
  static swig_type_info * swig_ECM_ptr  = SWIG_TypeQuery("Epetra_CrsMatrix       *");
  //static swig_type_info * swig_EMM_ptr  = SWIG_TypeQuery("Epetra_MsrMatrix       *");
  static swig_type_info * swig_EVM_ptr  = SWIG_TypeQuery("Epetra_VbrMatrix       *");
  static swig_type_info * swig_EVRM_ptr = SWIG_TypeQuery("Epetra_VbrRowMatrix    *");
  static swig_type_info * swig_EFVM_ptr = SWIG_TypeQuery("Epetra_FEVbrMatrix     *");
  static swig_type_info * swig_EFCM_ptr = SWIG_TypeQuery("Epetra_FECrsMatrix     *");
  static swig_type_info * swig_EJM_ptr  = SWIG_TypeQuery("Epetra_JadMatrix       *");
  //
  // Attempt to downcast to Epetra_VbrRowMatrix
  const Epetra_VbrRowMatrix * evrm = dynamic_cast< const Epetra_VbrRowMatrix* >(eo);
  if (evrm) return SWIG_NewPointerObj((void*)evrm, swig_EVRM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_FEVbrMatrix
  const Epetra_FEVbrMatrix * efvm = dynamic_cast< const Epetra_FEVbrMatrix* >(eo);
  if (efvm) return SWIG_NewPointerObj((void*)efvm, swig_EFVM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_FECrsMatrix
  const Epetra_FECrsMatrix * efcm = dynamic_cast< const Epetra_FECrsMatrix* >(eo);
  if (efcm) return SWIG_NewPointerObj((void*)efcm, swig_EFCM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_JadMatrix
  const Epetra_JadMatrix * ejm = dynamic_cast< const Epetra_JadMatrix* >(eo);
  if (ejm) return SWIG_NewPointerObj((void*)ejm, swig_EJM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_BasicRowMatrix
  const Epetra_BasicRowMatrix * ebrm = dynamic_cast< const Epetra_BasicRowMatrix* >(eo);
  if (ebrm) return SWIG_NewPointerObj((void*)ebrm, swig_EBRM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_CrsMatrix
  const Epetra_CrsMatrix * ecm = dynamic_cast< const Epetra_CrsMatrix* >(eo);
  if (ecm) return SWIG_NewPointerObj((void*)ecm, swig_ECM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_MsrMatrix
  // const Epetra_MsrMatrix * emm = dynamic_cast< const Epetra_MsrMatrix* >(eo);
  // if (emm) return SWIG_NewPointerObj((void*)emm, swig_EMM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_VbrMatrix
  const Epetra_VbrMatrix * evm = dynamic_cast< const Epetra_VbrMatrix* >(eo);
  if (evm) return SWIG_NewPointerObj((void*)evm, swig_EVM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_RowMatrix
  const Epetra_RowMatrix * erm = dynamic_cast< const Epetra_RowMatrix* >(eo);
  if (erm) return SWIG_NewPointerObj((void*)erm, swig_ERM_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_InvOperator
  const Epetra_InvOperator * eio = dynamic_cast< const Epetra_InvOperator* >(eo);
  if (eio) return SWIG_NewPointerObj((void*)eio, swig_EIO_ptr, cnvt_flags);
  //
  // Attempt to downcast to Epetra_FastCrsOperator
  // const Epetra_FastCrsOperator * efco = dynamic_cast< const Epetra_FastCrsOperator* >(eo);
  // if (efco) return SWIG_NewPointerObj((void*)efco, swig_EFCO_ptr, cnvt_flags);
  //
  // If downcast not possible, return python  Epetra_Operator
  return SWIG_NewPointerObj((void*)eo, swig_EO_ptr, cnvt_flags);
}

#endif

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS

Teuchos::RCP< const Epetra_Map >
PyTrilinos::getEpetraMapPtrFromEpetraBlockMap(const Epetra_BlockMap & ebm)
{
  const Epetra_Map * em_ptr  = dynamic_cast< const Epetra_Map* >(&ebm);
  if (!em_ptr)
  {
    PyErr_SetString(PyExc_TypeError, "Cannot upcast BlockMap to Map");
    throw PythonException();
  }
  return Teuchos::rcp(em_ptr, false);
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_Vector >
PyTrilinos::getEpetraVectorObjectAttr(PyObject * object, CONST char * name)
{
  static swig_type_info * swig_EV_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.Vector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP< Epetra_Vector > result = 
    *reinterpret_cast< Teuchos::RCP< Epetra_Vector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< Epetra_Vector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Epetra_Vector >
PyTrilinos::getConstEpetraVectorObjectAttr(PyObject * object, CONST char * name)
{
  static swig_type_info * swig_EV_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.Vector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP< const Epetra_Vector > result = 
    *reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Epetra_Vector >
PyTrilinos::getConstEpetraVectorItemObjectAttr(PyObject * object, CONST char * name, int i)
{
  static swig_type_info * swig_EV_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
  void * argp;
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(item, &argp, swig_EV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not tuple of type Epetra.Vector", name);
    Py_DECREF(item);
    throw PythonException();
  }
  Teuchos::RCP< const Epetra_Vector > result = 
    *reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  Py_DECREF(item);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_MultiVector >
PyTrilinos::getEpetraMultiVectorObjectAttr(PyObject * object, CONST char * name)
{
  static swig_type_info * swig_EMV_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_MultiVector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EMV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.MultiVector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_MultiVector > result = 
    *reinterpret_cast< Teuchos::RCP< Epetra_MultiVector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< Epetra_MultiVector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Epetra_MultiVector >
PyTrilinos::getConstEpetraMultiVectorObjectAttr(PyObject * object, CONST char * name)
{
  static swig_type_info * swig_EMV_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_MultiVector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EMV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.MultiVector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP< const Epetra_MultiVector > result = 
    *reinterpret_cast< Teuchos::RCP< const Epetra_MultiVector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< const Epetra_MultiVector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_Operator >
PyTrilinos::getEpetraOperatorObjectAttr(PyObject * object, CONST char * name)
{
  static swig_type_info * swig_EO_ptr = SWIG_TypeQuery("Teuchos::RCP< Epetra_Operator > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EO_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.Operator", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_Operator > result = 
    *reinterpret_cast< Teuchos::RCP< Epetra_Operator > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< Epetra_Operator > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

#endif   // HAVE_TEUCHOS
