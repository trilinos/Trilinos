// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Amesos_config.h"
#include "Amesos_Factory.h"
#include "Amesos_Klu.h"
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#include "Epetra_Object.h"


Amesos_BaseSolver* Amesos_Factory::Create( AmesosClassType ClassType, 
			     const Epetra_LinearProblem& LinearProblem, 
			     const Teuchos::ParameterList &ParameterList ) {

  switch( ClassType ) {
  case AMESOS_MUMPS:
#ifdef HAVE_AMESOS_MUMPS
    return new Amesos_Mumps(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Mumps is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  case AMESOS_SCALAPACK:
#ifdef HAVE_AMESOS_SCALAPACK
    return new Amesos_Scalapack(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Scalapack is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  case AMESOS_UMFPACK:
#ifdef HAVE_AMESOS_UMFPACK
    return new Amesos_Umfpack(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Umfpack is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  case AMESOS_DSCPACK:
#ifdef HAVE_AMESOS_DSCPACK
    return new Amesos_Dscpack(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Dscpack is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  case AMESOS_KLU:
#ifdef HAVE_AMESOS_KLU
    return new Amesos_Klu(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Klu is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  case AMESOS_SUPERLUDIST:
#ifdef HAVE_AMESOS_SUPERLUDIST
    return new Amesos_Superludist(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Superludist is not implemented" << endl ; 
    return 0 ; 
#endif
  case AMESOS_SUPERLU:
#ifdef HAVE_AMESOS_SUPERLU
    return new Amesos_Superlu(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Superlu is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  default:
    cerr << "Unknown class type" << endl ; 
    return 0 ; 
  }
}


