// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPFACTOROBJECT_HPP
#define IFPACK2_OVERLAPFACTOROBJECT_HPP

//! Ifpack2_OverlapFactorObject: Supports functionality common to Ifpack2 overlap factorization classes.

class Ifpack2_OverlapFactorObject {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Ifpack2_OverlapGraph.
  /*! Creates an object from the overlap graph. 
    \param In
           OverlapGraph - Graph describing the graph that should be used for the factors.
  */
  Ifpack2_OverlapFactorObject(const Ifpack2_OverlapGraph * OverlapGraph);

  //! Constructor using Tpetra_RowMatrix.
  /*! Creates an Ifpack2_Graph object from the user graph implicitly defined by the
	 Tpetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Tpetra_RowMatrix interface.
  */
  Ifpack2_OverlapFactorObject(const Tpetra_RowMatrix * UserMatrix);
  
  //! Copy constructor.
  Ifpack2_OverlapFactorObject(const Ifpack2_OverlapFactorObject & Source);

  //! Ifpack2_OverlapFactorObject Destructor
  virtual ~Ifpack2_OverlapFactorObject();
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
      be implemented by any derived class of Ifpack2_OverlapFactorObject.
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
  Ifpack2_OverlapGraph * OverlapGraph_;
  Tpetra_RowMatrix * UserMatrix_;
};
#endif // IFPACK2_OVERLAPFACTOROBJECT_HPP
