
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "NOX_Epetra_Interface.H"
#include "NOX_Utils.H"

#include "NOX_Epetra_FiniteDifferenceColoring.H"

using namespace NOX;
using namespace NOX::Epetra;

FiniteDifferenceColoring::FiniteDifferenceColoring(Interface& i, 
                             const Epetra_Vector& x, 
                             Epetra_MapColoring& colorMap_,
                             vector<Epetra_IntVector>& columns_,
                             double beta_, double alpha_) :
  FiniteDifference(i, x, beta_, alpha_),
  colorMap(&colorMap_),
  columns(&columns_)
{
  label = "NOX::FiniteDifferenceColoring Jacobian";
}

FiniteDifferenceColoring::FiniteDifferenceColoring(Interface& i, 
                             const Epetra_Vector& x, 
                             Epetra_CrsGraph& rawGraph_,
                             Epetra_MapColoring& colorMap_,
                             vector<Epetra_IntVector>& columns_,
                             double beta_, double alpha_) :
  FiniteDifference(i, x, rawGraph_, beta_, alpha_),
  colorMap(&colorMap_),
  columns(&columns_)
{
  label = "NOX::FiniteDifferenceColoring Jacobian";
}

FiniteDifferenceColoring::~FiniteDifferenceColoring()
{}

bool FiniteDifferenceColoring::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  // First check to make sure Jac is a NOX::Epetra::FiniteDifferenceColoring object
  FiniteDifferenceColoring* testMatrix = dynamic_cast<FiniteDifferenceColoring*>(&Jac);
  if (testMatrix == 0) {
    cout << "ERROR: NOX::Epetra::FiniteDifferenceColoring::computeJacobian() - "
	 << "Jacobian to evaluate is not a FiniteDifferenceColoring object!" << endl;
    throw "NOX Error";
  } 

  // We need the Epetra_CrsMatrix inside the FiniteDifferenceColoring object 
  // for the correct insertion commands.
  Epetra_CrsMatrix& jac = *testMatrix->jacobian;

  // Zero out Jacobian
  jacobian->PutScalar(0.0);

  int numColors = colorMap->NumColors();
  int* colorList = colorMap->ListOfColors();

  double eta = 0.0;  // Value to perturb the solution vector 
  eta = beta;
  
  int min = map.MinAllGID();  // Minimum Global ID value
  int max = map.MaxAllGID();  // Maximum Global ID value
  int myMin = map.MinMyGID(); // Minimum Local ID value
  int myMax = map.MaxMyGID(); // Maximum Local ID value

  // Compute the RHS at the initial solution
  interface.computeF(x, fo, Interface::Jacobian);

  x_perturb = x;    

  // loop over each color in the colorGraph
  for (int k = 0; k < numColors; k++) {

    // Perturb the solution vector using coloring
    Epetra_Map* cMap = colorMap->GenerateMap(colorList[k]);
    Epetra_Import Importer(map, *cMap);
    Epetra_Vector colorVect(*cMap);
    colorVect.PutScalar(eta);
    Epetra_Vector mappedColorVect(map);
    mappedColorVect.PutScalar(0.0);
    mappedColorVect.Import(colorVect, Importer, Insert);
    x_perturb.Update(1.0, mappedColorVect, 1.0);
    delete cMap; // clean up 

    // Compute the perturbed RHS
    interface.computeF(x_perturb, fp, Interface::Jacobian);
    
    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    Jc.Scale(1.0/eta);
   
    // Insert nonzero column entries into the jacobian    
    for (int j = myMin; j < myMax+1; j++) {
      int globalColumnID = columns->operator[](k)[j];
      if( globalColumnID >= 0) { 
	int err = jac.ReplaceGlobalValues(j,1,&Jc[map.LID(j)],&globalColumnID);
        if(err) {
          cout << "ERROR: Inserting global value with indices (" << j << ","
               << globalColumnID << ") = " << Jc[map.LID(j)] << endl;
        }
      }
    }

    // Unperturb the solution vector
    x_perturb = x;    

  }

  jac.TransformToLocal();

  return true;
}
