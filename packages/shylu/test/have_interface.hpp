
#ifndef HAVE_INTERFACE_HPP
#define HAVE_INTERFACE_HPP


#include <assert.h>
#include <iostream>
#include <sstream>


#include "Isorropia_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "Ifpack_ConfigDefs.h"
#include "shylu.h"
#include "shylu_util.h"
#include "Ifpack_ShyLU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

//namespace shylu_interface{

  using std::cout;
  using std::endl;
  
   using Teuchos::ParameterList;




template <class MatrixType>
int have_interface (MatrixType A, ParameterList shyLUList);
  
template <class scale>
int test(scale a);



#endif
//}
