//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef EDT_LINEARPROBLEM_ZOLTAN_H
#define EDT_LINEARPROBLEM_ZOLTAN_H

#include <Epetra_Transform.h>

class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_Export;
class Epetra_Import;

namespace EpetraExt {

///
/** Generates an Epetra_LinearProblem based on the repartitioning of an
 * original Epetra_LinearProblem using Zoltan partitioning algorithms
 */

class LinearProblem_Zoltan : public SameTypeTransform<Epetra_LinearProblem>
{
  const bool verbose_;

  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

  Epetra_LinearProblem * OldProblem_;
  Epetra_CrsGraph * OldGraph_;
  Epetra_CrsMatrix * OldMatrix_;
  Epetra_MultiVector * OldLHS_;
  Epetra_MultiVector * OldRHS_;
  Epetra_Map * OldRowMap_;

  Epetra_LinearProblem * NewProblem_;
  Epetra_CrsMatrix * NewMatrix_;
  Epetra_MultiVector * NewLHS_;
  Epetra_MultiVector * NewRHS_;

 public:

  ///
  /** Destructor
   */
  ~LinearProblem_Zoltan();

  ///
  /** Constructor
   */
  LinearProblem_Zoltan( bool verbose = false )
  : verbose_(verbose),
    OldProblem_(0),
    OldGraph_(0),
    OldMatrix_(0),
    OldLHS_(0),
    OldRHS_(0),
    OldRowMap_(0),
    NewProblem_(0),
    NewMatrix_(0),
    NewLHS_(0),
    NewRHS_(0),
    Importer_(0),
    Exporter_(0)
  {}

  ///
  /** Generates repartitioned Epetra_LinearProblem from original object
   */
  NewTypeRef operator()( OriginalTypeRef orig );

  ///
  /** Migrates data from original to new object.
   */
  bool fwd();

  ///
  /** Migrates data from new to original object.
   */
  bool rvs();

};

} //namespace EpetraExt

#endif //EDT_LINEARPROBLEM_ZOLTAN_H
