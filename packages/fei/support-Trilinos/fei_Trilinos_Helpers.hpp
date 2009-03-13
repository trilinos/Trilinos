/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Trilinos_Helpers_hpp_
#define _fei_Trilinos_Helpers_hpp_

#include "fei_trilinos_macros.hpp"
#include "fei_fwd.hpp"

#include <fei_Include_Trilinos.hpp>

#include <fei_mpi.h>
#include <fei_SharedPtr.hpp>

#include <fei_LinearProblemManager.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <fei_MatrixGraph.hpp>

namespace Trilinos_Helpers {

#ifdef HAVE_FEI_EPETRA

  /** Epetra_Map objects are light-weight wrt copying, since they employ a
      memory-model that uses a reference-counted pointer to an internal data
      object. Thus, it is reasonable for this function to return a pass-by-value
      Epetra_Map object.
  */
  Epetra_Map create_Epetra_Map(MPI_Comm comm,
                               const std::vector<int>& local_eqns);

  /** Epetra_BlockMap objects are light-weight wrt copying, since they employ a
      memory-model that uses a reference-counted pointer to an internal data
      object. Thus, it is reasonable for this function to return a pass-by-value
      Epetra_BlockMap object.
  */
  Epetra_BlockMap
    create_Epetra_BlockMap(const fei::SharedPtr<fei::VectorSpace>& vecspace);

  Epetra_CrsGraph
    create_Epetra_CrsGraph(const fei::SharedPtr<fei::MatrixGraph>& matgraph,
                           bool blockEntries,
                           bool orderRowsWithLocalColsFirst=false);

  fei::SharedPtr<fei::Matrix>
    create_from_Epetra_Matrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                              bool blockEntryMatrix,
                              fei::SharedPtr<fei::Reducer> reducer,
                              bool orderRowsWithLocalColsFirst=false);

  fei::SharedPtr<fei::Matrix>
    create_from_LPM_EpetraBasic(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                                 bool blockEntryMatrix,
                                 fei::SharedPtr<fei::Reducer> reducer,
                                 fei::SharedPtr<fei::LinearProblemManager>
                                   lpm_epetrabasic);
#endif

  /** Copies parameters from fei::ParameterSet to Teuchos::ParameterList.
    Does not clear any pre-existing contents from the Teuchos:ParameterList.
  */
  void copy_parameterset(const fei::ParameterSet& paramset,
                         Teuchos::ParameterList& paramlist);

  /** Copies parameters from Teuchos::ParameterList to fei::ParameterSet.
    Does not clear any pre-existing contents from the fei:ParameterSet.
  */
  void copy_parameterlist(const Teuchos::ParameterList& paramlist,
                          fei::ParameterSet& paramset);

#ifdef HAVE_FEI_EPETRA
  /** Extracts a pointer to a Epetra_MultiVector from a fei::Vector. Throws
    an exception if unsuccessful.
  */
  Epetra_MultiVector*
    get_Epetra_MultiVector(fei::Vector* feivec, bool soln_vec);

  /** Extracts a pointer to a Epetra_VbrMatrix from a fei::Matrix. Throws
    an exception if unsuccessful.
  */
  Epetra_VbrMatrix* get_Epetra_VbrMatrix(fei::Matrix* feimat);

  /** Extracts a pointer to a Epetra_CrsMatrix from a fei::Matrix. Throws
    an exception if unsuccessful.
  */
  Epetra_CrsMatrix* get_Epetra_CrsMatrix(fei::Matrix* feimat);

  /** Extracts pointers to epetra objects from fei container-objects.
    If epetra objects can't be obtained, output arguments are set
    to NULL.
  */
  void get_Epetra_pointers(fei::SharedPtr<fei::Matrix> feiA,
                           fei::SharedPtr<fei::Vector> feix,
                           fei::SharedPtr<fei::Vector> feib,
                           Epetra_CrsMatrix*& crsA,
                           Epetra_Operator*& opA,
                           Epetra_MultiVector*& x,
                           Epetra_MultiVector*& b);
#endif // HAVE_FEI_EPETRA

}//namespace Trilinos_Helpers

#endif // _Trilinos_Helpers_hpp_

