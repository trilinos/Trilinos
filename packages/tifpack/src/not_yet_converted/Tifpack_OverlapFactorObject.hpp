/*@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
*/

#ifndef TIFPACK_OVERLAPFACTOROBJECT_HPP
#define TIFPACK_OVERLAPFACTOROBJECT_HPP

//! Tifpack_OverlapFactorObject: Supports functionality common to Tifpack overlap factorization classes.

class Tifpack_OverlapFactorObject {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Tifpack_OverlapGraph.
  /*! Creates an object from the overlap graph. 
    \param In
           OverlapGraph - Graph describing the graph that should be used for the factors.
  */
  Tifpack_OverlapFactorObject(const Tifpack_OverlapGraph * OverlapGraph);

  //! Constructor using Tpetra_RowMatrix.
  /*! Creates an Tifpack_Graph object from the user graph implicitly defined by the
	 Tpetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Tpetra_RowMatrix interface.
  */
  Tifpack_OverlapFactorObject(const Tpetra_RowMatrix * UserMatrix);
  
  //! Copy constructor.
  Tifpack_OverlapFactorObject(const Tifpack_OverlapFactorObject & Source);

  //! Tifpack_OverlapFactorObject Destructor
  virtual ~Tifpack_OverlapFactorObject();
  //@}

  //@{ \name Initialization methods.

  //! Initialize values from user matrix A, can be called repeatedly as matrix values change.
  /*! Processes matrix values, primarily handling overlap if any has been requested.  This method
      then calls ProcessOverlapMatrix(), a virtual method that must be implemented by any class
      that derives from this class.
    \param In 
           UserMatrix - User matrix to be processed.
   */
  virtual int InitValues(const Tpetra_RowMatrix * UserMatrix);

  //! Compute factors.
  /*! This function computes factors using the method DerivedFactor() that 
      is implemented by the derived class.
    InitValues() must be called before the factorization can proceed.
   */
  virtual int Factor();
  //@}

  //@{ \name Attribue accessor methods.


  //! If storage has been allocated, this query returns true, otherwise it returns false.
  bool Allocated() const {return(Allocated_);};

  //! If values have been initialized, this query returns true, otherwise it returns false.
  bool ValuesInitialized() const {return(ValuesInitialized_);};

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool Factored() const {return(Factored_);};
   //@}
 
 protected:

  //@{ \name Methods that must be implemented by derived classes.
  //! Virtual method that processes the overlap matrix as needed by the derived class.
  /*! This method is called by InitValues() afer the user matrix has been distributed 
      to support overlap (if any overlap is requested).  ProcessOverlapMatrix must
      be implemented by any derived class of Tifpack_OverlapFactorObject.
  */
  virtual int ProcessOverlapMatrix(const Tpetra_RowMatrix &A)=0;

  //! Virtual method that computes the factors as needed by the derived class.
  /*! This method is called by Factor() afer some safety checks have been performed.
  */
  virtual int DerivedFactor()=0;
   //@}

  void SetAllocated(bool Flag) {Allocated_ = Flag;};
  void SetFactored(bool Flag) {Factored_ = Flag;};
  void SetValuesInitialized(bool Flag) {ValuesInitialized_ = Flag;};

  bool Factored_;
  bool Allocated_;
  bool ValuesInitialized_;
  Tifpack_OverlapGraph * OverlapGraph_;
  Tpetra_RowMatrix * UserMatrix_;
};
#endif // TIFPACK_OVERLAPFACTOROBJECT_HPP
