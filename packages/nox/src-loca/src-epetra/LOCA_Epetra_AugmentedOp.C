// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_config.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"

#include "LOCA_Epetra_AugmentedOp.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Epetra::AugmentedOp::AugmentedOp(
 const Teuchos::RCP<LOCA::GlobalData>& global_data,
 const Teuchos::RCP<Epetra_Operator>& jac,
 const Teuchos::RCP<const Epetra_MultiVector>& a_,
 const Teuchos::RCP<const Epetra_MultiVector>& b_,
 const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> c_)
  : globalData(global_data),
    label("LOCA::Epetra::AugmentedOp"),
    jacOperator(jac),
    underlyingMap(jacOperator->OperatorDomainMap()),
    underlyingComm(underlyingMap.Comm()),
    localMap(a_->NumVectors(), 0, jacOperator->Comm()),
    a(a_),
    b(b_),
    c(View, localMap, c_->values(), c_->stride(), c_->numCols()),
    underlyingLength(a->MyLength()),
    numConstraints(a->NumVectors()),
    extendedMapPtr(NULL),
    extendedImportMapPtr(NULL),
    extendedImporter(NULL),
    importedInput(NULL),
    result_y(NULL),
    tmp(NULL),
    haveParamComponent(false),
    useTranspose(false),
    dlapack()
{

  // Determine if this processor stores parameter component
  haveParamComponent = (underlyingComm.MyPID() == 0);

  // Build extended map
  buildExtendedMap(underlyingMap, extendedMapPtr, false,
           haveParamComponent);

  // Build extended importer map
  buildExtendedMap(underlyingMap, extendedImportMapPtr, true,
           haveParamComponent);

  // Build importer
  extendedImporter = new Epetra_Import(*extendedImportMapPtr,
                       *extendedMapPtr);
}

LOCA::Epetra::AugmentedOp::~AugmentedOp()
{
  delete extendedMapPtr;
  delete extendedImportMapPtr;
  delete extendedImporter;
  delete importedInput;
  delete result_y;
  delete tmp;
}

int
LOCA::Epetra::AugmentedOp::SetUseTranspose(bool UseTranspose)
{
  useTranspose = UseTranspose;
  jacOperator->SetUseTranspose(UseTranspose);
  return 0;
}

int
LOCA::Epetra::AugmentedOp::Apply(const Epetra_MultiVector& Input,
                 Epetra_MultiVector& Result) const
{

  // Get number of vectors
  int n = Input.NumVectors();

  // Check num vectors is correct
  if (importedInput == NULL || n != importedInput->NumVectors())
    globalData->locaErrorCheck->throwError("LOCA::Epetra::AugmentedOp::Apply()"
                       "Must call init() before Apply()!");

  // Import parameter components
  importedInput->Import(Input, *extendedImporter, Insert);

  // get views
  double **input_view;
  double **result_view;
  importedInput->ExtractView(&input_view);
  Result.ExtractView(&result_view);

  // get views of paramter components
  // this is done by setting the pointer for each column to start at
  // underlyingLength
  double **input_view_y = new double*[n];
  for (int i=0; i<n; i++)
    input_view_y[i] = input_view[i]+underlyingLength;

  // break input, result into components
  Epetra_MultiVector input_x(View, underlyingMap, input_view, n);
  Epetra_MultiVector input_y(View, localMap, input_view_y, n);
  Epetra_MultiVector result_x(View, underlyingMap, result_view, n);

  // Compute J*input_x
  jacOperator->Apply(input_x, result_x);

  if (useTranspose) {

    // Compute J^T*input_x + b*input_y
    result_x.Multiply('N', 'N', 1.0, *b, input_y, 1.0);

    // Compute a^T*input_x + c^T*input_y
    result_y->Multiply('T', 'N', 1.0, *a, input_x, 0.0);
    result_y->Multiply('T', 'N', 1.0, c, input_y, 1.0);

  }
  else {

    // Compute J*input_x + a*input_y
    result_x.Multiply('N', 'N', 1.0, *a, input_y, 1.0);

    // Compute b^T*input_x + c*input_y
    result_y->Multiply('T', 'N', 1.0, *b, input_x, 0.0);
    result_y->Multiply('N', 'N', 1.0, c, input_y, 1.0);

  }

  // Set parameter component
  if (haveParamComponent)
    for (int j=0; j<n; j++)
      for (int i=0; i<numConstraints; i++)
    result_view[j][underlyingLength+i] = (*result_y)[j][i];

   delete [] input_view_y;

  return 0;
}

int
LOCA::Epetra::AugmentedOp::ApplyInverse(const Epetra_MultiVector& Input,
                    Epetra_MultiVector& Result) const
{

  // Get number of vectors
  int n = Input.NumVectors();

  // Check num vectors is correct
  if (importedInput == NULL || n != importedInput->NumVectors())
    globalData->locaErrorCheck->throwError(
                 "LOCA::Epetra::AugmentedOp::ApplyInverse()"
                 "Must call init() before ApplyInverse()!");

  // Import parameter components
  importedInput->Import(Input, *extendedImporter, Insert);

  // get views
  double **input_view;
  double **result_view;
  importedInput->ExtractView(&input_view);
  Result.ExtractView(&result_view);

  // get views of paramter components
  // this is done by setting the pointer for each column to start at
  // underlyingLength
  double **input_view_y = new double*[n];
  for (int i=0; i<n; i++)
    input_view_y[i] = input_view[i]+underlyingLength;

  // break input, result into components
  Epetra_MultiVector input_x(View, underlyingMap, input_view, n);
  Epetra_MultiVector input_y(View, localMap, input_view_y, n);
  Epetra_MultiVector result_x(View, underlyingMap, result_view, n);

  // copy input_y into result_y
  result_y->Scale(1.0, input_y);

  // Temporary epetra vectors
  Epetra_MultiVector tmp2(c);
  tmp->PutScalar(0.0);

  // Solve J*result_x = input_x
  jacOperator->ApplyInverse(input_x, result_x);

  if (useTranspose) {

    // Solve J*tmp = b
    jacOperator->ApplyInverse(*b, *tmp);

    // Compute input_y - a^T*result_x
    result_y->Multiply('T', 'N', -1.0, *a, result_x, 1.0);

    // Compute c - a^T*tmp
    tmp2.Multiply('T', 'N', -1.0, *a, *tmp, 1.0);

  }
  else {

    // Solve J*tmp = a
    jacOperator->ApplyInverse(*a, *tmp);

    // Compute input_y - b^T*result_x
    result_y->Multiply('T', 'N', -1.0, *b, result_x, 1.0);

    // Compute c - b^T*tmp1
    tmp2.Multiply('T', 'N', -1.0, *b, *tmp, 1.0);

  }

  // Solve (c - z^T*tmp1)^-1 * (input_y - z^T*result_x)  where z = a or b
  int *ipiv = new int[numConstraints];
  int info;
  double *result_y_view;
  int result_y_lda;
  result_y->ExtractView(&result_y_view, &result_y_lda);
  double *tmp2_view;
  int tmp2_lda;
  tmp2.ExtractView(&tmp2_view, &tmp2_lda);
  dlapack.GESV(numConstraints, n, tmp2_view, tmp2_lda, ipiv, result_y_view,
           result_y_lda, &info);
  delete [] ipiv;
  if (info != 0) {
    globalData->locaErrorCheck->throwError(
                 "LOCA::Epetra::AugmentedOp::ApplyInverse()"
                 "Solve of dense matrix failed!");
  }

  // Compute result_x = result_x - tmp*result_y
  result_x.Multiply('N', 'N', -1.0, *tmp, *result_y, 1.0);

  // Set parameter component
  if (haveParamComponent)
    for (int j=0; j<n; j++)
      for (int i=0; i<numConstraints; i++)
    result_view[j][underlyingLength+i] = (*result_y)[j][i];

  delete [] input_view_y;

  return 0;
}

double
LOCA::Epetra::AugmentedOp::NormInf() const
{
  double Jn, an, bn;
  double *t = new double[numConstraints];

  Jn = jacOperator->NormInf();

  a->NormInf(t);
  an = t[0];
  for (int i=1; i<numConstraints; i++)
    if (t[i] > an)
      an = t[i];

  b->NormInf(t);
  bn = t[0];
  for (int i=1; i<numConstraints; i++)
    if (t[i] > bn)
      bn = t[i];

  delete [] t;

  return Jn + an + bn;
}


const char*
LOCA::Epetra::AugmentedOp::Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
LOCA::Epetra::AugmentedOp::UseTranspose() const
{
  return useTranspose;
}

bool
LOCA::Epetra::AugmentedOp::HasNormInf() const
{
  return jacOperator->HasNormInf();
}

const Epetra_Comm&
LOCA::Epetra::AugmentedOp::Comm() const
{
  return underlyingComm;
}
const Epetra_Map&
LOCA::Epetra::AugmentedOp::OperatorDomainMap() const
{
  return *extendedMapPtr;
}

const Epetra_Map&
LOCA::Epetra::AugmentedOp::OperatorRangeMap() const
{
  return *extendedMapPtr;
}

void
LOCA::Epetra::AugmentedOp::init(const Epetra_MultiVector& x)
{
  if (importedInput == NULL || importedInput->NumVectors() != x.NumVectors()) {
    if (importedInput != NULL) {
      delete importedInput;
      delete result_y;
      delete tmp;
    }
    importedInput = new Epetra_MultiVector(*extendedImportMapPtr,
                       x.NumVectors());
    result_y = new Epetra_MultiVector(localMap, x.NumVectors());
    tmp = new Epetra_MultiVector(underlyingMap, x.NumVectors());
  }
}

Teuchos::RCP<Epetra_MultiVector>
LOCA::Epetra::AugmentedOp::buildEpetraAugmentedMultiVec(
                 const Epetra_MultiVector& x,
                 const NOX::Abstract::MultiVector::DenseMatrix *y,
                 bool doCopy) const
{
  Teuchos::RCP<Epetra_MultiVector> extVec =
    Teuchos::rcp(new Epetra_MultiVector(*extendedMapPtr, x.NumVectors()));

  if (doCopy) {
    for (int j=0; j<x.NumVectors(); j++) {
      for (int i=0; i<underlyingLength; i++)
    (*extVec)[j][i] = x[j][i];

      if (haveParamComponent && y != NULL)
    for (int i=0; i<numConstraints; i++)
      (*extVec)[j][underlyingLength+i] = (*y)(i,j);
    }
  }
  return extVec;
}

void
LOCA::Epetra::AugmentedOp::setEpetraAugmentedMultiVec(
                 Epetra_MultiVector& x,
                 NOX::Abstract::MultiVector::DenseMatrix& y,
                 const Epetra_MultiVector& augMultiVec) const
{
  // Check num vectors is correct
  if (importedInput == NULL || x.NumVectors() != importedInput->NumVectors())
    globalData->locaErrorCheck->throwError(
             "LOCA::Epetra::AugmentedOp::setEpetraAugmentedMultiVec()"
             "Must call init() before setEpetraAugmentedMultiVec()!");

  // Import parameter components
  importedInput->Import(augMultiVec, *extendedImporter, Insert);

  for (int j=0; j<x.NumVectors(); j++) {
    for (int i=0; i<underlyingLength; i++)
      x[j][i] = (*importedInput)[j][i];
    for (int i=0; i<numConstraints; i++)
      y(i,j) = (*importedInput)[j][underlyingLength+i];
  }
}

void
LOCA::Epetra::AugmentedOp::buildExtendedMap(const Epetra_BlockMap& uMap,
                        Epetra_Map*& eMapPtr,
                        bool buildImporter,
                        bool haveParam)
{
  Epetra_BlockMap& nonconstUnderlyingMap = const_cast<Epetra_BlockMap&>(uMap);

  // Convert underlying map to point map if necessary
  Epetra_Map* uPointMapPtr =
    dynamic_cast<Epetra_Map*>(&nonconstUnderlyingMap);
  bool allocatedPointMap = false;
  if (uPointMapPtr == NULL) {
    allocatedPointMap = true;
    blockMap2PointMap(uMap, uPointMapPtr);
  }

  int max_gid = uPointMapPtr->MaxAllGID();
  int num_global_elements = uPointMapPtr->NumGlobalElements();
  int num_my_elements = uPointMapPtr->NumMyElements();
  int *global_elements = uPointMapPtr->MyGlobalElements();
  const Epetra_Comm& comm = uPointMapPtr->Comm();
  int index_base = uPointMapPtr->IndexBase();

  int ext_num_global_elements;
  int ext_num_my_elements;
  int *ext_global_elements;

  // Compute number of extended global elements
  if (buildImporter)
    ext_num_global_elements =
      num_global_elements + numConstraints*comm.NumProc();
  else
    ext_num_global_elements = num_global_elements + numConstraints;

  // Compute number of extended local elements
  if (buildImporter || haveParam)
     ext_num_my_elements = num_my_elements + numConstraints;
  else
    ext_num_my_elements = num_my_elements;

  // Allocate extended global elements array
  ext_global_elements = new int[ext_num_my_elements];

  // Set extended global elements
  for (int i=0; i<num_my_elements; i++) {
    ext_global_elements[i] = global_elements[i];
  }
  if (buildImporter || haveParam)
    for (int i=0; i<numConstraints; i++)
      ext_global_elements[num_my_elements+i] = max_gid + 1 + i;

  // Create extended point map
  eMapPtr = new Epetra_Map(ext_num_global_elements, ext_num_my_elements,
               ext_global_elements, index_base, comm);

  // Free global elements array
  delete [] ext_global_elements;
  if (allocatedPointMap)
    delete uPointMapPtr;
}

int
LOCA::Epetra::AugmentedOp::blockMap2PointMap(const Epetra_BlockMap& BlockMap,
                    Epetra_Map*& PointMap) const
{
  // Generate an Epetra_Map that has the same number and distribution of points
  // as the input Epetra_BlockMap object.  The global IDs for the output PointMap
  // are computed by using the MaxElementSize of the BlockMap.  For variable block
  // sizes this will create gaps in the GID space, but that is OK for Epetra_Maps.

  int MaxElementSize = BlockMap.MaxElementSize();
  int PtNumMyElements = BlockMap.NumMyPoints();

  TEUCHOS_ASSERT_INEQUALITY(PtNumMyElements, >, 0);
  int * PtMyGlobalElements = new int[PtNumMyElements];

  int NumMyElements = BlockMap.NumMyElements();

  int curID = 0;
  for (int i=0; i<NumMyElements; i++) {
    int StartID = BlockMap.GID(i)*MaxElementSize;
    int ElementSize = BlockMap.ElementSize(i);
    for (int j=0; j<ElementSize; j++) PtMyGlobalElements[curID++] = StartID+j;
  }
  assert(curID==PtNumMyElements); // Sanity test

  PointMap = new Epetra_Map(-1, PtNumMyElements, PtMyGlobalElements, BlockMap.IndexBase(), BlockMap.Comm());

  if (PtNumMyElements>0) delete [] PtMyGlobalElements;

  // Maps not compatible?
  TEUCHOS_TEST_FOR_EXCEPT(!BlockMap.PointSameAs(*PointMap));

  return(0);
}
