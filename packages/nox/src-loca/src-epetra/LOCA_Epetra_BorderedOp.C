//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "Epetra_config.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"

#include "NOX_Utils.H"
#include "LOCA_Epetra_BorderedOp.H"

LOCA::Epetra::BorderedOp::BorderedOp(Epetra_Operator& jac, 
				     const Epetra_Vector& a, 
				     const Epetra_Vector& b) :
  label("LOCA::Epetra::BorderedOp"),
  jacOperator(jac),
  aVecPtr(&a),
  bVecPtr(&b),
  underlyingMapPtr(NULL),
  underlyingCommPtr(NULL),
  underlyingLength(0),
  extendedMapPtr(NULL),
  extendedImportMapPtr(NULL),
  extendedImporter(NULL),
  haveParamComponent(false),
  useTranspose(false)
{
  // Get underlying map
  //underlyingMapPtr = &(a.Map());
  underlyingMapPtr = &(jac.OperatorDomainMap());

  // Get underlying communicator
  underlyingCommPtr = &(underlyingMapPtr->Comm());

  // Get underlying local vector length
  underlyingLength = a.MyLength();

  // Determine if this processor stores parameter component
  haveParamComponent = (underlyingCommPtr->MyPID() == 0);

  // Build extended map
  buildExtendedMap(*underlyingMapPtr, extendedMapPtr, false, 
		   haveParamComponent);

  // Build extended importer map
  buildExtendedMap(*underlyingMapPtr, extendedImportMapPtr, true, 
		   haveParamComponent);

  // Build importer
  extendedImporter = new Epetra_Import(*extendedImportMapPtr,
				       *extendedMapPtr);
}

LOCA::Epetra::BorderedOp::~BorderedOp()
{
  delete extendedMapPtr;
  delete extendedImportMapPtr;
  delete extendedImporter;
}

int 
LOCA::Epetra::BorderedOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  jacOperator.SetUseTranspose(UseTranspose);
  return 0;
}

int 
LOCA::Epetra::BorderedOp::Apply(const Epetra_MultiVector& cInput, 
			Epetra_MultiVector& Result) const
{
  // cast-away const on cInput
  Epetra_MultiVector& Input = const_cast<Epetra_MultiVector&>(cInput);

  // cast input and result to vectors
  Epetra_Vector& input = dynamic_cast<Epetra_Vector&>(Input);
  Epetra_Vector& result = dynamic_cast<Epetra_Vector&>(Result);

  // Create imported input vector
  Epetra_Vector input_import(*extendedImportMapPtr);

  // Import paramter component
  input_import.Import(input, *extendedImporter, Insert);

  // get views
  double *input_view;
  double *result_view;
  input_import.ExtractView(&input_view);
  result.ExtractView(&result_view);

  // break input, result into components
  Epetra_Vector input_x(View, *underlyingMapPtr, input_view);
  double input_y = input_view[underlyingLength];

  Epetra_Vector result_x(View, *underlyingMapPtr, result_view);
  double w;

  // Compute J*input_x
  jacOperator.Apply(input_x, result_x);

  if (useTranspose) {
    
    // Compute J*input_x + input_y*b
    result_x.Update(input_y, *bVecPtr, 1.0);

    // Compute a^T*input_x
    aVecPtr->Dot(input_x, &w);

  }
  else {
    
    // Compute J*input_x + input_y*a
    result_x.Update(input_y, *aVecPtr, 1.0);

    // Compute b^T*input_x
    bVecPtr->Dot(input_x, &w);

  }

  // Set parameter component
  if (haveParamComponent)
    result_view[underlyingLength] = w;

  // Note:  we don't need to explicitly set the solution (x) and 
  // null vector (y) components of result since the epetra vectors 
  // were created with views.

  return 0;
}

int 
LOCA::Epetra::BorderedOp::ApplyInverse(const Epetra_MultiVector& cInput, 
			       Epetra_MultiVector& Result) const
{
  // cast-away const on cInput
  Epetra_MultiVector& Input = const_cast<Epetra_MultiVector&>(cInput);

  // cast input and result to vectors
  Epetra_Vector& input = dynamic_cast<Epetra_Vector&>(Input);
  Epetra_Vector& result = dynamic_cast<Epetra_Vector&>(Result);

  // Create imported input vector
  Epetra_Vector input_import(*extendedImportMapPtr);

  // Import paramter component
  input_import.Import(input, *extendedImporter, Insert);

  // get views
  double *input_view;
  double *result_view;
  input_import.ExtractView(&input_view);
  result.ExtractView(&result_view);

   // break input, result into components
  Epetra_Vector input_x(View, *underlyingMapPtr, input_view);
  double input_y = input_view[underlyingLength];

  Epetra_Vector result_x(View, *underlyingMapPtr, result_view);
  double u, v, w;

  // Temporary epetra vectors
  Epetra_Vector tmp(*underlyingMapPtr); tmp.PutScalar(0.0);
 
  // Solve J*result_x = input_x 
  jacOperator.ApplyInverse(input_x, result_x);

  if (useTranspose) {

    // Solve J*tmp = b
    jacOperator.ApplyInverse(*bVecPtr, tmp);

    // Compute a^T*result_x
    aVecPtr->Dot(result_x, &u);

    // Compute a^T*tmp
    aVecPtr->Dot(tmp, &v);

  }
  else {

    // Solve J*tmp = a
    jacOperator.ApplyInverse(*aVecPtr, tmp);

    // Compute b^T*result_x
    bVecPtr->Dot(result_x, &u);

    // Compute b^T*tmp
    bVecPtr->Dot(tmp, &v);
  }

  w = (u - input_y)/v;
  result_x.Update(-w, tmp, 1.0);

  if (haveParamComponent)
    result_view[underlyingLength] = w;

  return 0;
}

double 
LOCA::Epetra::BorderedOp::NormInf() const
{
  double Jn, an, bn;

  Jn = jacOperator.NormInf();
  aVecPtr->NormInf(&an);
  bVecPtr->NormInf(&bn);

  return Jn + an + bn;
}


char* 
LOCA::Epetra::BorderedOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
LOCA::Epetra::BorderedOp::UseTranspose() const
{
  return useTranspose;
}

bool 
LOCA::Epetra::BorderedOp::HasNormInf() const
{
  return jacOperator.HasNormInf();
}

const Epetra_Comm & 
LOCA::Epetra::BorderedOp::Comm() const
{
  return *underlyingCommPtr;
}
const Epetra_Map& 
LOCA::Epetra::BorderedOp::OperatorDomainMap() const
{
  return *extendedMapPtr;
}

const Epetra_Map& 
LOCA::Epetra::BorderedOp::OperatorRangeMap() const
{
  return *extendedMapPtr;
}

Epetra_Vector*
LOCA::Epetra::BorderedOp::buildEpetraExtendedVec(Epetra_Vector& x, double p,
					 bool doCopy) const
{
  Epetra_Vector* extVec = new Epetra_Vector(*extendedMapPtr);

  if (doCopy) {
    double *x_view; x.ExtractView(&x_view);
    double *ext_view; extVec->ExtractView(&ext_view);

    for (int i=0; i<underlyingLength; i++) {
      ext_view[i] = x_view[i];
    }
    if (haveParamComponent)
      ext_view[underlyingLength] = p;
  }

  return extVec;
}

void
LOCA::Epetra::BorderedOp::setEpetraExtendedVec(Epetra_Vector& x,double& p,
				       Epetra_Vector& extVec) const
{
  // Import parameter component
  Epetra_Vector extVecImport(*extendedImportMapPtr);
  extVecImport.Import(extVec, *extendedImporter, Insert);

  double *x_view; x.ExtractView(&x_view);
  double *ext_view; extVecImport.ExtractView(&ext_view);
  
  for (int i=0; i<underlyingLength; i++) {
    x_view[i] = ext_view[i];
  }
  p = ext_view[underlyingLength];
}

void
LOCA::Epetra::BorderedOp::buildExtendedMap(const Epetra_BlockMap& uMap,
				   Epetra_Map*& eMapPtr,
				   bool buildImporter,
				   bool haveParam)
{
  Epetra_BlockMap& nonconstUnderlyingMap = const_cast<Epetra_BlockMap&>(uMap);

  // Convert underlying map to point map if necessary
  Epetra_Map* uPointMapPtr = 
    dynamic_cast<Epetra_Map*>(&nonconstUnderlyingMap);
  if (uPointMapPtr == NULL)
    blockMap2PointMap(uMap, uPointMapPtr);

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
    ext_num_global_elements = num_global_elements + comm.NumProc();
  else
    ext_num_global_elements = num_global_elements + 1;

  // Compute number of extended local elements
  if (buildImporter || haveParam)
     ext_num_my_elements = num_my_elements + 1;
  else
    ext_num_my_elements = num_my_elements;

  // Allocate extended global elements array
  ext_global_elements = new int[ext_num_my_elements];

  // Set extended global elements
  for (int i=0; i<num_my_elements; i++) {
    ext_global_elements[i] = global_elements[i];
  }
  if (buildImporter || haveParam)
    ext_global_elements[num_my_elements] = max_gid + 1;

  // Create extended point map
  eMapPtr = new Epetra_Map(ext_num_global_elements, ext_num_my_elements,
			   ext_global_elements, index_base, comm);

  // Free global elements array
  delete [] ext_global_elements;
}

int
LOCA::Epetra::BorderedOp::blockMap2PointMap(const Epetra_BlockMap& BlockMap,
				    Epetra_Map*& PointMap) const
{
  // Generate an Epetra_Map that has the same number and distribution of points
  // as the input Epetra_BlockMap object.  The global IDs for the output PointMap
  // are computed by using the MaxElementSize of the BlockMap.  For variable block
  // sizes this will create gaps in the GID space, but that is OK for Epetra_Maps.

  int MaxElementSize = BlockMap.MaxElementSize();
  int PtNumMyElements = BlockMap.NumMyPoints();
  int * PtMyGlobalElements = 0;
  if (PtNumMyElements>0) PtMyGlobalElements = new int[PtNumMyElements];

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

  if (!BlockMap.PointSameAs(*PointMap)) {EPETRA_CHK_ERR(-1);} // Maps not compatible
  return(0);
}
