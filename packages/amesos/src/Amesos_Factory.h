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

#ifndef _AMESOS_FACTORY_H_
#define _AMESOS_FACTORY_H_
#include "Amesos_BaseSolver.h"
#include "AmesosClassType.h"

//! Amesos_Factory:  A method for creating Amesos classes
/*!  Amesos_Factory allows a code to delay the decision about which
concrete class to use to implement the Amesos_BaseSolver interface.  
*/
//
class Amesos_Factory { 

  //@{ \name Creation method
  //! Amesos_Factory Create method
  /*! Creates an instance of the Amesos_BaseSolver class specified by 
    ClassType 
  */
public: 
  Amesos_BaseSolver *Create( AmesosClassType ClassType, 
			     const Epetra_LinearProblem& LinearProblem, 
			     const AMESOS::Parameter::List &ParameterList );
};  // End of  class Amesos_Factory  
#endif /* _AMESOS_FACTORY_H_ */
