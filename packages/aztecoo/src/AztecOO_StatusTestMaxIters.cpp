
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
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

#include "AztecOO_StatusTestMaxIters.h"

AztecOO_StatusTestMaxIters::AztecOO_StatusTestMaxIters(int MaxIters) 
  : AztecOO_StatusTest() {

  if (MaxIters < 1) MaxIters_ = 1;
    
  MaxIters_ = MaxIters;
  status_ = Unchecked;
}

AztecOO_StatusType AztecOO_StatusTestMaxIters::CheckStatus(int CurrentIter, 
							   Epetra_MultiVector * CurrentResVector, 
							   double CurrentResNormEst,
							   bool SolutionUpdated) {
  status_ = Unconverged;
  Niters_ = CurrentIter;
  if (Niters_ >= MaxIters_)
    status_ = Failed;
  return status_;
}


ostream& AztecOO_StatusTestMaxIters::Print(ostream& stream, int indent) const {

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  PrintStatus(stream, status_);
  stream << "Number of Iterations = ";
  stream << Niters_;
  stream << ((Niters_<MaxIters_) ? " < " : ((Niters_==MaxIters_) ? " = " : " > "));
  stream << MaxIters_;
  stream << endl;
 return stream;
}
