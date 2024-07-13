// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xfer_Operator.H"
#include "GenericEpetraProblem.H"
#include "Epetra_Vector.h"

#undef DEBUG_TRANSFER_OPERATOR

// Constructor
XferOp::XferOp(GenericEpetraProblem& probA, const GenericEpetraProblem& probB)
{
  // Establish mappings from problem meshes.  This is based solely on each
  // problem's mesh and involves a simple linear interpolation (possibly
  // extrapolation).
  Epetra_Vector &meshA = probA.getMesh();
  Epetra_Vector &meshB = const_cast<GenericEpetraProblem&>(probB).getMesh();

  // For now, this is designed for serial execution since some problemA
  // nodes will lie within off-processor problemB elements in general.

  // The following assumes that 1D elements exist such that an element is
  // defined by a left and right end-node which are sequentially numbered.
  double tol(1.e-8);
  for( int i = 0; i < meshA.MyLength(); i++)
  {
    double xlocA = meshA[i];
    for( int j = 0; j < (meshB.MyLength()-1); ++j )
    {
      double xlocB_left  = meshB[j]   ;
      double xlocB_right = meshB[j+1] ;
      // First check for what we choose to define as a direct match (within tol)
      if( fabs(xlocB_left - xlocA) <= tol )
      {
        dependentNodes.insert( multimap<int,int>::value_type(i, j) );
        dependentWeights.insert( multimap<int,double>::value_type(i, double(1.0)) );
        break;
      }
      // Next determine which two nodes of meshB xlocA lies between
      else if( (xlocA > xlocB_left) && (xlocA < xlocB_right) )
      {
        dependentNodes.insert( multimap<int,int>::value_type(i, j) );
        dependentNodes.insert( multimap<int,int>::value_type(i, j+1) );
        double wt_right = (xlocA - xlocB_left)/(xlocB_right - xlocB_left);
        double wt_left = 1.0 - wt_right;
        dependentWeights.insert( multimap<int,double>::value_type(i, wt_left) );
        dependentWeights.insert( multimap<int,double>::value_type(i, wt_right) );
        break;
      }
      // Finally check for a direct match of the right node
      else if( fabs(xlocB_right - xlocA) <= tol )
      {
        dependentNodes.insert( multimap<int,int>::value_type(i, j+1) );
        dependentWeights.insert( multimap<int,double>::value_type(i, double(1.0)) );
        break;
      }
    }
    if( !dependentNodes.count(i) )
      std::cout << "WARNING: Node " << i << " at position " << xlocA
           << " not interpolated to transfer mesh !!" << std::endl;
  }

#ifdef DEBUG_TRANSFER_OPERATOR
  for( int i = 0; i < meshA.MyLength(); i++)
  {
    std::cout << "meshA node " << i << " at loc --> " << meshA[i]
         << " needs node(s) from meshB: ";
    pair< multimap<int, int>::iterator,
          multimap<int, int>::iterator > rangeN
          = dependentNodes.equal_range(i);
    pair< multimap<int, double>::iterator,
          multimap<int, double>::iterator > rangeW
          = dependentWeights.equal_range(i);
    multimap<int, int>::iterator iterN;
    multimap<int, double>::iterator iterW;
    int j;
    for( j = 0, iterN = rangeN.first, iterW = rangeW.first;
           iterN != rangeN.second; j++, iterN++, iterW++)
      std::cout << "(" << (*iterN).second << ", " << (*iterW).second << ")   ";
    std::cout << std::endl;
  }
#endif

}

//-----------------------------------------------------------------------------

// Destructor
XferOp::~XferOp()
{
}

//-----------------------------------------------------------------------------

void
XferOp::transferField(Epetra_Vector& vecTo, Epetra_Vector& vecFrom)
{
  // Do the transfer using Epetra_Vectors
  for( int i = 0; i < vecTo.MyLength(); i++)
  {
    vecTo[i] = 0.0;
    pair< multimap<int, int>::iterator,
          multimap<int, int>::iterator > rangeN
          = dependentNodes.equal_range(i);
    pair< multimap<int, double>::iterator,
          multimap<int, double>::iterator > rangeW
          = dependentWeights.equal_range(i);
    multimap<int, int>::iterator iterN;
    multimap<int, double>::iterator iterW;
    int j;
    for( j = 0, iterN = rangeN.first, iterW = rangeW.first;
           iterN != rangeN.second; j++, iterN++, iterW++)
      vecTo[i] += (*iterW).second * vecFrom[(*iterN).second];
  }
}

//-----------------------------------------------------------------------------
