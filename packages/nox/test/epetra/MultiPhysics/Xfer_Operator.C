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
  for( int i = 0; i < meshA.MyLength(); i++) {
    double xlocA = meshA[i];
    for( int j = 0; j < (meshB.MyLength()-1); j++) {
      double xlocB_left = meshB[j];
      double xlocB_right = meshB[j+1];
      // First check for what we choose to define as a direct match (within tol)
      if( fabs(xlocB_left - xlocA) <= tol ) {
        dependentNodes.insert( pair<int,int>(i, j) );
        dependentWeights.insert( pair<int,double>(i, double(1.0)) );
        break;
      }
      // Next determine which two nodes of meshB xlocA lies between
      else if( (xlocA > xlocB_left) && (xlocA < xlocB_right) ) {
        dependentNodes.insert( pair<int,int>(i, j) );
        dependentNodes.insert( pair<int,int>(i, j+1) );
        double wt_right = (xlocA - xlocB_left)/(xlocB_right - xlocB_left);
        double wt_left = 1.0 - wt_right;
        dependentWeights.insert( pair<int,double>(i, wt_left) );
        dependentWeights.insert( pair<int,double>(i, wt_right) );
        break;
      }
      // Finally check for a direct match of the right node
      else if( fabs(xlocB_right - xlocA) <= tol ) {
        dependentNodes.insert( pair<int,int>(i, j+1) );
        dependentWeights.insert( pair<int,double>(i, double(1.0)) );
        break;
      }
    }
    if( !dependentNodes.count(i) )
      cout << "ERROR: Node " << i << " at position " << xlocA 
           << " not interpolated to transfer mesh !!" << endl;
  }

#ifdef DEBUG_TRANSFER_OPERATOR
  for( int i = 0; i < meshA.MyLength(); i++) {
    cout << "meshA node " << i << " at loc --> " << meshA[i]
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
      cout << "(" << iterN->second << ", " << iterW->second << ")   ";
    cout << endl;
  }
#endif
    
}

// Destructor
XferOp::~XferOp() { };

// Calculates the values of u and x at the specified gauss point
void XferOp::transferField(Epetra_Vector& vecTo, Epetra_Vector& vecFrom)
{
  // Do the transfer using Epetra_Vectors
  for( int i = 0; i < vecTo.MyLength(); i++) {
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
      vecTo[i] += iterW->second * vecFrom[iterN->second];
  }
}
