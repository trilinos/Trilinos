
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
