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
                                                                                
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_Utils.H"

#include "NOX_EpetraNew_FiniteDifferenceColoring.H"

using namespace NOX;
using namespace NOX::EpetraNew;

// This constructor is needed for inheritance but is inadequate for using
// coloring in parallel since the raw matrix graph is not known.
FiniteDifferenceColoring::FiniteDifferenceColoring(Interface::Required& i, 
                             const Epetra_Vector& x, 
                             Epetra_MapColoring& colorMap_,
                             vector<Epetra_IntVector>& columns_,
			     bool parallelColoring,
                             double beta_, double alpha_) :
  FiniteDifference(i, x, beta_, alpha_),
  coloringType(SERIAL),
  colorMap(&colorMap_),
  columns(&columns_),
  numColors(colorMap->NumColors()),
  maxNumColors(colorMap->MaxNumColors()),
  colorList(colorMap->ListOfColors()),
  cMap(0),
  Importer(0),
  colorVect(0),
  betaColorVect(0),
  mappedColorVect(0),
  columnMap(0),
  rowColImporter(0),
  xCol_perturb(0)
{
  label = "NOX::FiniteDifferenceColoring Jacobian";

  if( parallelColoring )
    coloringType = PARALLEL;

  createColorContainers();
}

FiniteDifferenceColoring::FiniteDifferenceColoring(Interface::Required& i, 
                             const Epetra_Vector& x, 
                             Epetra_CrsGraph& rawGraph_,
                             Epetra_MapColoring& colorMap_,
                             vector<Epetra_IntVector>& columns_,
			     bool parallelColoring,
                             double beta_, double alpha_) :
  FiniteDifference(i, x, rawGraph_, beta_, alpha_),
  coloringType(SERIAL),
  colorMap(&colorMap_),
  columns(&columns_),
  numColors(colorMap->NumColors()),
  maxNumColors(colorMap->MaxNumColors()),
  colorList(colorMap->ListOfColors()),
  cMap(0),
  Importer(0),
  colorVect(0),
  betaColorVect(0),
  columnMap(&rawGraph_.ColMap()),
  rowColImporter(new Epetra_Import(*columnMap, map)),
  xCol_perturb(new Epetra_Vector(rawGraph_.ColMap())),
  mappedColorVect(new Epetra_Vector(rawGraph_.ColMap()))
{
  label = "NOX::FiniteDifferenceColoring Jacobian";

  if( parallelColoring ) 
    coloringType = PARALLEL;

  createColorContainers();
}

FiniteDifferenceColoring::~FiniteDifferenceColoring()
{
  delete cMap; cMap = 0;
  delete Importer; Importer = 0;
  delete colorVect; colorVect = 0;
  delete betaColorVect; betaColorVect = 0;
  delete rowColImporter; rowColImporter = 0;
  delete xCol_perturb; xCol_perturb = 0;
  delete mappedColorVect; mappedColorVect = 0;
}

bool FiniteDifferenceColoring::computeJacobian(const Epetra_Vector& x)
{
  return( computeJacobian(x, *this));
}

bool FiniteDifferenceColoring::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  // First check to make sure Jac is a 
  // NOX::EpetraNew::FiniteDifferenceColoring object
  FiniteDifferenceColoring* testMatrix = 
         dynamic_cast<FiniteDifferenceColoring*>(&Jac);
  if (testMatrix == 0) {
    cout << "ERROR: NOX::EpetraNew::FiniteDifferenceColoring::computeJacobian() - "
	 << "Jacobian to evaluate is not a FiniteDifferenceColoring object!" 
         << endl;
    throw "NOX Error";
  } 

  // We need the Epetra_CrsMatrix inside the FiniteDifferenceColoring object 
  // for the correct insertion commands.
  Epetra_CrsMatrix& jac = *testMatrix->jacobian;

  // Zero out Jacobian
  jac.PutScalar(0.0);

  // Create an extra perturbed residual vector pointer if needed
  if ( diffType == Centered )
    if ( !fmPtr )
      fmPtr = new Epetra_Vector(x);

  // Create a reference to the extra perturbed residual vector
  Epetra_Vector& fm = *fmPtr;

  double scaleFactor = 1.0;
  if ( diffType == Backward )
    scaleFactor = -1.0;

  int min = map.MinAllGID();  // Minimum Global ID value
  int max = map.MaxAllGID();  // Maximum Global ID value
  int myMin = map.MinMyGID(); // Minimum Local ID value
  int myMax = map.MaxMyGID(); // Maximum Local ID value

  // We need to loop over the largest number of colors on a processor
  
  // Use the overlap (column-space) version of the solution
  xCol_perturb->Import(x, *rowColImporter, Insert);
  
  // Compute the RHS at the initial solution
  if( coloringType == SERIAL )
    computeF(*xCol_perturb, fo, NOX::EpetraNew::Interface::Required::FD_Res);
  else {
//    x_perturb.Export(*xCol_perturb, *rowColImporter, Insert);
    for (int i=0; i<x_perturb.MyLength(); i++)
      x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID(i))];
    computeF(x_perturb, fo, NOX::EpetraNew::Interface::Required::FD_Res);
  }
  
  // loop over each color in the colorGraph
  list<int>::iterator allBegin = listOfAllColors.begin(),
                      allEnd   = listOfAllColors.end(),
                      allIter;
  
  std::map<int, int>::iterator mapEnd = colorToNumMap.end(),
                               myMapIter;
  for ( allIter = allBegin; allIter != allEnd; ++allIter ) {
    bool skipIt = true;       // Used to screen colors this proc does not have
    int k = -1;               // index in colorList of color
    myMapIter = colorToNumMap.find( *allIter );
    if( myMapIter != mapEnd ) {
      skipIt = false;
      k = myMapIter->second;
    }

    // Perturb the solution vector using coloring
  
    // ----- First create color map and assoc perturbation vectors
    int color = -1;
    if( !skipIt ) color = colorList[k];
    cMap = colorMap->GenerateMap(color);
    colorVect = new Epetra_Vector(*cMap);
    betaColorVect = new Epetra_Vector(*cMap);
    betaColorVect->PutScalar(beta);
   
    // ----- Fill colorVect with computed perturbation values
    // NOTE that we do the mapping ourselves here instead of using an
    // Epetra_Import object.  This is to ensure local mapping only since
    // Import/Export operations involve off-processor transfers, which we 
    // wish to avoid.
    for (int i=0; i<colorVect->MyLength(); i++)
      (*colorVect)[i] = (*xCol_perturb)[columnMap->LID(cMap->GID(i))];
    colorVect->Abs(*colorVect);
    colorVect->Update(1.0, *betaColorVect, alpha);

    // ----- Map perturbation vector to original index space 
    // ----- and use it 
    mappedColorVect->PutScalar(0.0);
    // Here again we do the mapping ourselves to avoid off-processor data
    // transfers that would accompany use of Epetra_Import/Export objects.
    for (int i=0; i<colorVect->MyLength(); i++)
      (*mappedColorVect)[columnMap->LID(cMap->GID(i))] = (*colorVect)[i];
    xCol_perturb->Update(scaleFactor, *mappedColorVect, 1.0);

    // Compute the perturbed RHS
    if( coloringType == SERIAL )
      computeF(*xCol_perturb, fp, NOX::EpetraNew::Interface::Required::FD_Res);
    else {
      //x_perturb.Export(*xCol_perturb, *rowColImporter, Insert);
      for (int i=0; i<x_perturb.MyLength(); i++)
        x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID(i))];
      computeF(x_perturb, fp, NOX::EpetraNew::Interface::Required::FD_Res);
    }
    
    if ( diffType == Centered ) {
      xCol_perturb->Update(-2.0, *mappedColorVect, 1.0);
      if( coloringType == SERIAL )
        computeF(*xCol_perturb, fm, 
                 NOX::EpetraNew::Interface::Required::FD_Res);
      else {
        //x_perturb.Export(*xCol_perturb, *rowColImporter, Insert);
        for (int i=0; i<x_perturb.MyLength(); i++)
          x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID(i))];
        computeF(x_perturb, fm, NOX::EpetraNew::Interface::Required::FD_Res);
      }
    }

    // Compute the column k of the Jacobian
    // Note that division by the perturbation is delayed until insertion below
    if ( diffType != Centered ) {
      Jc.Update(1.0, fp, -1.0, fo, 0.0);
    }
    else {
      Jc.Update(1.0, fp, -1.0, fm, 0.0);
    }
    
    // Insert nonzero column entries into the jacobian    
    if( !skipIt ) {
      for (int j = myMin; j < myMax+1; j++) {
        int globalColumnID = (*columns)[k][map.LID(j)];
        if( globalColumnID >= 0) { 
  
          // Now complete the approximation to the derivative by dividing by
          // the appropriate perturbation
          if ( diffType != Centered ) 
            Jc[map.LID(j)] /= 
              (scaleFactor * (*mappedColorVect)[columnMap->LID(globalColumnID)]);
          else
            Jc[map.LID(j)] /= 
              (2.0 * (*mappedColorVect)[columnMap->LID(globalColumnID)]);
  
  	int err = jac.ReplaceGlobalValues(j,1,&Jc[map.LID(j)],&globalColumnID);
          if(err) {
            cout << "ERROR (" << map.Comm().MyPID() << ") : "
                 << "Inserting global value with indices (" << j << ","
                 << globalColumnID << ") = " << Jc[map.LID(j)] << endl;
          }
        }
      }
    }

    // Clean up memory for color-dependent objects
    delete Importer; Importer = 0;
    delete betaColorVect; betaColorVect = 0;
    delete colorVect; colorVect = 0;
    delete cMap; cMap = 0;

    // Unperturb the solution vector
    xCol_perturb->Import(x, *rowColImporter, Insert);
  }

  // Use a barrier to be safe
  x.Comm().Barrier();

  jac.TransformToLocal();

//  jac.Print(cout);

  return true;
}

void FiniteDifferenceColoring::createColorContainers()
{
  // First send all procs' color ids to all others
  int sumNumColors = numColors;
  colorMap->Comm().SumAll(&numColors, &sumNumColors, 1);
  int * allColorList = new int[maxNumColors*colorMap->Comm().NumProc()];
  int * myColorList = new int[maxNumColors];
  for( int i=0; i< maxNumColors; i++)
    if( i<numColors )
      myColorList[i] = colorList[i];
    else
      myColorList[i] = -1;

  colorMap->Comm().GatherAll( myColorList, allColorList, maxNumColors );
  
  // Insert all colors into our list
  for( int i = 0; i < maxNumColors*colorMap->Comm().NumProc(); i++ )
    listOfAllColors.push_back( allColorList[i] );
  
  listOfAllColors.remove(-1);
  listOfAllColors.sort();
  listOfAllColors.unique();
  
  // Now create a map to use with the index object
  for( int i = 0; i < numColors; i++ )
    colorToNumMap.insert( pair< int, int > ( colorList[i], i ) );
  
}
