/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_OVERLAPFACTOROBJECT_H
#define IFPACK_OVERLAPFACTOROBJECT_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

//! Ifpack_OverlapFactorObject: Supports functionality common to Ifpack overlap factorization classes.

class Ifpack_OverlapFactorObject {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Ifpack_OverlapGraph.
  /*! Creates an object from the overlap graph. 
    \param In
           OverlapGraph - Graph describing the graph that should be used for the factors.
  */
  Ifpack_OverlapFactorObject(const Ifpack_OverlapGraph * OverlapGraph);

  //! Constructor using Epetra_RowMatrix.
  /*! Creates an Ifpack_Graph object from the user graph implicitly defined by the
	 Epetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Epetra_RowMatrix interface.
  */
  Ifpack_OverlapFactorObject(const Epetra_RowMatrix * UserMatrix);
  
  //! Copy constructor.
  Ifpack_OverlapFactorObject(const Ifpack_OverlapFactorObject & Source);

  //! Ifpack_OverlapFactorObject Destructor
  virtual ~Ifpack_OverlapFactorObject();
  //@}

  //@{ \name Initialization methods.

  //! Initialize values from user matrix A, can be called repeatedly as matrix values change.
  /*! Processes matrix values, primarily handling overlap if any has been requested.  This method
      then calls ProcessOverlapMatrix(), a virtual method that must be implemented by any class
      that derives from this class.
    \param In 
           UserMatrix - User matrix to be processed.
   */
  virtual int InitValues(const Epetra_RowMatrix * UserMatrix);

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
      be implemented by any derived class of Ifpack_OverlapFactorObject.
  */
  virtual int ProcessOverlapMatrix(const Epetra_RowMatrix &A)=0;

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
  Ifpack_OverlapGraph * OverlapGraph_;
  Epetra_RowMatrix * UserMatrix_;
};
#endif // IFPACK_OVERLAPFACTOROBJECT_H
