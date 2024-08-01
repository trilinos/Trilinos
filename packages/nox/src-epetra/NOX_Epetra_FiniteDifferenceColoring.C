// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
#include "Epetra_IntVector.h"
#else
#include "Epetra_LongLongVector.h"
#endif
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "Epetra_Time.h" // For timing performance

#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Utils.H"

#include "NOX_Epetra_FiniteDifferenceColoring.H"

using namespace NOX;
using namespace NOX::Epetra;

// This constructor is needed for inheritance but is inadequate for using
// coloring in parallel since the raw matrix graph is not known.
FiniteDifferenceColoring::FiniteDifferenceColoring(
     Teuchos::ParameterList& printingParams,
     const Teuchos::RCP<Interface::Required>& i,
     const NOX::Epetra::Vector& x,
     const Teuchos::RCP<Epetra_MapColoring>& colorMap_,
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
     const Teuchos::RCP<std::vector<Epetra_IntVector> >& columns_,
#else
     const Teuchos::RCP<std::vector<Epetra_LongLongVector> >& columns_,
#endif
     bool parallelColoring,
     bool distance1_,
     double beta_, double alpha_) :
  FiniteDifference(printingParams, i, x, beta_, alpha_),
  coloringType(NOX_SERIAL),
  distance1(distance1_),
  colorMap(colorMap_),
  columns(columns_),
  numColors(colorMap->NumColors()),
  maxNumColors(colorMap->MaxNumColors()),
  colorList(colorMap->ListOfColors()),
  cMap(0),
  Importer(0),
  colorVect(0),
  betaColorVect(0),
  mappedColorVect(0),
  xCol_perturb(0),
  columnMap(0),
  rowColImporter(0)
{
  label = "NOX::FiniteDifferenceColoring Jacobian";

  if( parallelColoring )
    coloringType = NOX_PARALLEL;

  createColorContainers();
}

FiniteDifferenceColoring::FiniteDifferenceColoring(
         Teuchos::ParameterList& printingParams,
     const Teuchos::RCP<Interface::Required>& i,
     const NOX::Epetra::Vector& x,
     const Teuchos::RCP<Epetra_CrsGraph>& rawGraph_,
     const Teuchos::RCP<Epetra_MapColoring>& colorMap_,
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
     const Teuchos::RCP<std::vector<Epetra_IntVector> >& columns_,
#else
     const Teuchos::RCP<std::vector<Epetra_LongLongVector> >& columns_,
#endif
     bool parallelColoring,
     bool distance1_,
     double beta_, double alpha_) :
  FiniteDifference(printingParams, i, x, rawGraph_, beta_, alpha_),
  coloringType(NOX_SERIAL),
  distance1(distance1_),
  colorMap(colorMap_),
  columns(columns_),
  numColors(colorMap->NumColors()),
  maxNumColors(colorMap->MaxNumColors()),
  colorList(colorMap->ListOfColors()),
  cMap(0),
  Importer(0),
  colorVect(0),
  betaColorVect(0),
  mappedColorVect(new Epetra_Vector(rawGraph_->ColMap())),
  xCol_perturb(new Epetra_Vector(rawGraph_->ColMap())),
  columnMap(&(rawGraph_->ColMap())),
  rowColImporter(new Epetra_Import(*columnMap, fo.Map()))
{
  label = "NOX::FiniteDifferenceColoring Jacobian";

  if( parallelColoring )
    coloringType = NOX_PARALLEL;

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
  // NOX::Epetra::FiniteDifferenceColoring object
  FiniteDifferenceColoring* testMatrix =
         dynamic_cast<FiniteDifferenceColoring*>(&Jac);
  if (testMatrix == 0) {
    std::cout << "ERROR: NOX::Epetra::FiniteDifferenceColoring::computeJacobian() - "
     << "Jacobian to evaluate is not a FiniteDifferenceColoring object!"
         << std::endl;
    throw std::runtime_error("NOX Error");
  }

  const Epetra_BlockMap& map = fo.Map();

  // Create a timer for performance
  Epetra_Time fillTimer(x.Comm());

  // We need the Epetra_CrsMatrix inside the FiniteDifferenceColoring object
  // for the correct insertion commands.
  Epetra_CrsMatrix& jac = *testMatrix->jacobian;

  // Zero out Jacobian
  jac.PutScalar(0.0);

  // Create an extra perturbed residual vector pointer if needed
  if ( diffType == Centered )
    if ( Teuchos::is_null(fmPtr) )
      fmPtr = Teuchos::rcp(new Epetra_Vector(x));

  double scaleFactor = 1.0;
  if ( diffType == Backward )
    scaleFactor = -1.0;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int myMin = map.MinMyGID(); // Minimum locally owned GID
  int myMax = map.MaxMyGID(); // Maximum locally owned GID
#else
  long long myMin = map.MinMyGID64(); // Minimum locally owned GID
  long long myMax = map.MaxMyGID64(); // Maximum locally owned GID
#endif

  // We need to loop over the largest number of colors on a processor

  // Use the overlap (column-space) version of the solution
  xCol_perturb->Import(x, *rowColImporter, Insert);

  // Compute the RHS at the initial solution
  if( coloringType == NOX_SERIAL )
    computeF(*xCol_perturb, fo, NOX::Epetra::Interface::Required::FD_Res);
  else {
//    x_perturb.Export(*xCol_perturb, *rowColImporter, Insert);
    for (int i=0; i<x_perturb.MyLength(); i++)
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID(i))];
#else
      x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID64(i))];
#endif
    computeF(x_perturb, fo, NOX::Epetra::Interface::Required::FD_Res);
  }

  // loop over each color in the colorGraph
  std::list<int>::iterator allBegin = listOfAllColors.begin(),
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
      k = (*myMapIter).second;
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      (*colorVect)[i] = (*xCol_perturb)[columnMap->LID(cMap->GID(i))];
#else
      (*colorVect)[i] = (*xCol_perturb)[columnMap->LID(cMap->GID64(i))];
#endif
    colorVect->Abs(*colorVect);
    colorVect->Update(1.0, *betaColorVect, alpha);

    // ----- Map perturbation vector to original index space
    // ----- and use it
    mappedColorVect->PutScalar(0.0);
    // Here again we do the mapping ourselves to avoid off-processor data
    // transfers that would accompany use of Epetra_Import/Export objects.
    for (int i=0; i<colorVect->MyLength(); i++)
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      (*mappedColorVect)[columnMap->LID(cMap->GID(i))] = (*colorVect)[i];
#else
      (*mappedColorVect)[columnMap->LID(cMap->GID64(i))] = (*colorVect)[i];
#endif
    xCol_perturb->Update(scaleFactor, *mappedColorVect, 1.0);

    // Compute the perturbed RHS
    if( coloringType == NOX_SERIAL )
      computeF(*xCol_perturb, fp, NOX::Epetra::Interface::Required::FD_Res);
    else {
      //x_perturb.Export(*xCol_perturb, *rowColImporter, Insert);
      for (int i=0; i<x_perturb.MyLength(); i++)
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
        x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID(i))];
#else
        x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID64(i))];
#endif
      computeF(x_perturb, fp, NOX::Epetra::Interface::Required::FD_Res);
    }

    if ( diffType == Centered ) {
      xCol_perturb->Update(-2.0, *mappedColorVect, 1.0);
      if( coloringType == NOX_SERIAL )
        computeF(*xCol_perturb, *fmPtr,
                 NOX::Epetra::Interface::Required::FD_Res);
      else {
        //x_perturb.Export(*xCol_perturb, *rowColImporter, Insert);
        for (int i=0; i<x_perturb.MyLength(); i++)
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
          x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID(i))];
#else
          x_perturb[i] = (*xCol_perturb)[columnMap->LID(map.GID64(i))];
#endif
        computeF(x_perturb, *fmPtr, NOX::Epetra::Interface::Required::FD_Res);
      }
    }

    // Compute the column k of the Jacobian
    // Note that division by the perturbation is delayed until insertion below
    if ( diffType != Centered ) {
      Jc.Update(1.0, fp, -1.0, fo, 0.0);
    }
    else {
      Jc.Update(1.0, fp, -1.0, *fmPtr, 0.0);
    }

    // Insert nonzero column entries into the jacobian
    if( !skipIt ) {
      for (int j = myMin; j < myMax+1; j++) {
        // Allow for the possibility that rows j from myMin to myMax are not necessarily contigous
        if (!map.MyGID(j))
          continue;
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
        int globalColumnID = (*columns)[k][map.LID(j)];
#else
        long long globalColumnID = (*columns)[k][map.LID(j)];
#endif

        // If using distance1 coloring, only allow diagonal fills
        if( distance1 && (j != globalColumnID) )
          continue; // Skip off-diagonals

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
            std::cout << "ERROR (" << map.Comm().MyPID() << ") : "
                 << "Inserting global value with indices (" << j << ","
                 << globalColumnID << ") = " << Jc[map.LID(j)] << std::endl;
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

  double fillTime = fillTimer.ElapsedTime();
  x.Comm().Barrier();

  if (utils.isPrintType(Utils::Details)) {
    for(int n = 0; n < map.Comm().NumProc(); n++) {
      if(map.Comm().MyPID() == n)
        std::cout << "\tTime to fill Jacobian [" << n << "] --> "
             << fillTime << " sec." << std::endl;
      x.Comm().Barrier();
    }
    std::cout << std::endl;
  }

  jac.FillComplete();

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
    colorToNumMap[ colorList[i] ] =  i;

  // Cleanup
  delete [] myColorList;
  delete [] allColorList;
}
