
/* Copyright (2003) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2003, the United States Government is granted for itself and
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

#include "Amesos_config.h"
#include "Amesos_Factory.h"
#include "Amesos_Klu.h"
#include "Epetra_Object.h"

Amesos_BaseSolver* Amesos_Factory::Create( AmesosClassType ClassType, 
			     const Epetra_LinearProblem& LinearProblem, 
			     const AMESOS::Parameter::List &ParameterList ) {

  switch( ClassType ) {
  case AMESOS_KLU:
#ifdef HAVE_AMESOS_KLU
    return new Amesos_Klu(LinearProblem,ParameterList); 
#else
    cerr << "Amesos_Klu is not implemented" << endl ; 
    return 0 ; 
#endif
    break;
  default:
    cerr << "Unknown class type" << endl ; 
    return 0 ; 
  }
}
