
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

#ifndef _EPETRA_SRCDISTOBJECT_H_
#define _EPETRA_SRCDISTOBJECT_H_
class Epetra_BlockMap;


//! Epetra_SrcDistObject: A class for supporting flexible source distributed objects for import/export operations.

/*! The Epetra_SrcDistObject is a base class for all Epetra distributed global objects that are potential 
    source objects for the general Epetra_DistObject class.  It provides a way to send a very general distributed
    object as the potential source object for an import or export object.  For example, it is possible to pass
    an Epetra_RowMatrix object as the source object for an import/export where the target is an Epetra_CrsMatrix, or
    an Epetra_CrsGraph (where the RowMatrix values will be ignored).

*/

//==========================================================================
class Epetra_SrcDistObject {

  public:
  //@{ \name Destructor.
  //! Epetra_SrcDistObject destructor.  
  virtual ~Epetra_SrcDistObject() {};
  //@}

  
  //@{ \name Attribute accessor methods.
  //! Returns the address of the Epetra_BlockMap for this multi-vector.
  virtual const Epetra_BlockMap & Map() const = 0;
};

#endif /* _EPETRA_SRCDISTOBJECT_H_ */
