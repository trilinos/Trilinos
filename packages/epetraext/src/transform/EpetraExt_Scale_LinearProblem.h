//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
                                                                                                    
#ifndef EpetraExt_LINEARPROBLEM_SCALE_H
#define EpetraExt_LINEARPROBLEM_SCALE_H

#include <vector>

#include <EpetraExt_Transform.h>

class Epetra_LinearProblem;
class Epetra_Vector;

namespace EpetraExt {

///
/** Given an input Epetra_LinearProblem, recursive, left and right scaling are performed.
 */
class LinearProblem_Scale : public InPlaceTransform<Epetra_LinearProblem>
{
 public:

  enum ScaleType{ Sum, Max, Diag, None };

  ///
  /** Destructor
   */
  ~LinearProblem_Scale();

  ///
  /** Constructor
   */
  LinearProblem_Scale( ScaleType left = Sum,
                       ScaleType right = Sum,
                       double exp_fac = 1.0,
                       int iterations = 1 )
  : lScale_(left),
    rScale_(right),
    expFac_(exp_fac),
    iters_(iterations),
    scaled_(false)
  {}

  ///
  /** Applies forward scaling
   */
  bool fwd();

  ///
  /** Reverses scaling
   */
  bool rvs();

 private:

  const ScaleType lScale_;
  const ScaleType rScale_;

  const double expFac_;

  const int iters_;

  bool scaled_;

  std::vector<Epetra_Vector*> lScaleVecs_;
  std::vector<Epetra_Vector*> rScaleVecs_;
};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_SCALE_H

