/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
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

#ifndef IFPACK_OVERLAPFACTOROBJECT_H
#define IFPACK_OVERLAPFACTOROBJECT_H

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
  Ifpack_OverlapFactorObject(const Ifpack_CrsIlut & Source);

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

  //@{ \name Methods that must be implemented by derived classes
.
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

  bool Allocated() const {return(Allocated_);};

  bool Factored_;
  bool Allocated_;
  bool ValuesInitialized_;
  Ifpack_OverlapGraph * OverlapGraph_;
  Epetra_RowMatrix * UserMatrix_;
};
#endif // IFPACK_OVERLAPFACTOROBJECT_H
