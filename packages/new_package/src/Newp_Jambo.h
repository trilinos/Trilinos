
/* Copyright (2003) Sandia Corportation. Under the terms of Contract 
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

#ifndef _NEWP_JAMBO_H_
#define _NEWP_JAMBO_H_
#include "Epetra_ConfigDefs.h"
#include "Newp_Jambo.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
//! Newp_Jambo: A sample class 

/*! The Newp_Jambo class prints out a "Hello World" type message using the only word of 
Swahili that I know: Jambo, which means Hello.  The sole purpose of this code is to demonstrate
the conditional compilation and linkage features of the AUTOTOOLS package of tools.

<b>A typical heading</b>
<ul>
  <li> A typical first list entry
  <li> A typical second list entry
</ul>

<b>Another typical heading</b>

*/

//=========================================================================
class Newp_Jambo {

  public:

  //@{ \name Constructors/destructors.
  //! Basic Newp_Jambo constuctor.
  /*! Creates a Newp_Jambo object and fills with default values.  

    \param In 
           Comm - A Epetra Communicator 

	   \warning This is a typical warning.

    \return Pointer to a Newp_Jambo.

  */
  Newp_Jambo(const Epetra_Comm& Comm);

  //! Newp_Jambo copy constructor.
  
  Newp_Jambo(const Newp_Jambo& Source);
  
  //@}
  
  //@{ \name Print methods


  //! Print method
  virtual void Print(ostream & os) const;
  //@}


 private:

  const Epetra_Comm& Comm_ ; 

};

#endif /* _NEWP_JAMBO_H_ */
