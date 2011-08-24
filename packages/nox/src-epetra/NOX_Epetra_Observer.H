// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef NOX_EPETRA_OBSERVER
#define NOX_EPETRA_OBSERVER

#include "Epetra_Vector.h"


namespace NOX {
namespace Epetra {

class Observer 
{
public:
  Observer () {};
  virtual ~Observer () {};

  //! Method called by Piro NOXSolver as a hook for postprocessing
  virtual void observeSolution(const Epetra_Vector& soln) = 0;

  //! LOCA calls this version, which will discard param info in this default implementation
  virtual void observeSolution(const Epetra_Vector& soln, double time_or_param_val)
  { observeSolution(soln); };

private:

};

}
}

#endif //NOX_EPETRA_NOXOBSERVER
