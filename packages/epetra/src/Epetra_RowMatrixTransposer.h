
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

#ifndef EPETRA_CRSMATRIXTRANSPOSER_H
#define EPETRA_CRSMATRIXTRANSPOSER_H
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_Map;
class Epetra_Export;

//! Epetra_RowMatrixTransposer: A class for transposing an Epetra_RowMatrix object.

/*! This class provides capabilities to construct a transpose matrix of an existing Epetra_RowMatrix
	  object and (optionally) redistribute it across a parallel distributed memory machine.
*/

class Epetra_RowMatrixTransposer {
    
  public:

  //@{ \name Constructors/destructors.
  //! Primary Epetra_RowMatrixTransposer constructor.
  /*!
    \param Matrix (In) An existing Epetra_RowMatrix object.  The Epetra_RowMatrix, the LHS and RHS pointers
		       do not need to be defined before this constructor is called.

    \return Pointer to a Epetra_RowMatrixTransposer object.

  */ 
  Epetra_RowMatrixTransposer(Epetra_RowMatrix * OrigMatrix);

  //! Epetra_RowMatrixTransposer copy constructor.
  
  Epetra_RowMatrixTransposer(const Epetra_RowMatrixTransposer& Source);
  
  //! Epetra_RowMatrixTransposer destructor.
  
  virtual ~Epetra_RowMatrixTransposer();
  //@}
  
  //@{ \name Forward transformation methods.
  
  //! Generate a new Epetra_CrsMatrix as the transpose of an Epetra_RowMatrix passed into the constructor.
  /*! Constructs a new Epetra_CrsMatrix that is a copy of the Epetra_RowMatrix passed in to the constructor.
		
		\param MakeDataContiguous (In) Causes the output matrix, LHS and RHS to be stored in a form compatible with
		       Fortran-style solvers.  The output matrix will be compatible with the Harwell-Boeing compressed
					 column format.  The RHS and LHS will be stored such that the last value in column j of the 
					 multivector is stored next to the first value in column j+1.
		\param TransposeRowMap (Optional/In) If this argument is defined, the transpose matrix will be distributed
		       using this map as the row map for the transpose.  If it is set to zero, the transpose matrix will use
					 the OrigMatrix->RowMatrixDomainMap as the row map.

		\return Integer error code, 0 if no errors.  Negative if some fatal error occured.
					 
  */
  int CreateTranspose(const bool MakeDataContiguous,
											Epetra_CrsMatrix *& TransposeMatrix,
											Epetra_Map * TransposeRowMap = 0);

	
  //! Update the values of an already-redistributed problem.
  /*! Updates the values of an already-redistributed problem.  This method allows updating 
		  the redistributed problem without
		  allocating new storage.

    \param MatrixWithNewValues (In) The values from MatrixWithNewValues will be copied into the TransposeMatrix.  The
		       MatrixWithNewValues object must be identical in structure to the original matrix object used to create
					 this instance of Epetra_RowMatrixTransposer.

		\return Integer error code, 0 if no errors.  Negative if some fatal error occured.
					 
  */
  int UpdateTransposeValues(Epetra_RowMatrix * MatrixWithNewValues);
  //@}
  
  //@{ \name Reverse transformation methods.
  //! Update values of original matrix (Not implemented and not sure if we will implement this).
   int UpdateOriginalMatrixValues();
  //@}
  
  //@{ \name Attribute accessor methods.

  //! Returns const reference to the Epetra_Map object describing the row distribution of the transpose matrix.
  /*! The RedistExporter object can be used to redistribute other Epetra_DistObject objects whose maps are compatible with
		  the original linear problem map, or with the RedistMap().
			\warning Must not be called before CreateTranspose()is called.
  */
  const Epetra_Map & TransposeRowMap() const {return(*TransposeRowMap_);};
  //! Returns const reference to the Epetra_Export object used to redistribute the original matrix.
  /*! The TransposeExporter object can be used to redistribute other Epetra_DistObject objects whose maps are compatible with
		  the original matrix. 
			\warning Must not be called before CreateTranspose() is called.
  */
  const Epetra_Export & TransposeExporter() const{return(*TransposeExporter_);};
  //@}
  
 private: 
	void DeleteData();
	Epetra_RowMatrix * OrigMatrix_;
	Epetra_CrsMatrix * TransposeMatrix_;
	Epetra_Export * TransposeExporter_;
	Epetra_Map * TransposeRowMap_;
	bool MakeDataContiguous_;
	bool TransposeCreated_;
	int NumMyRows_;
	int NumMyCols_;
	int MaxNumEntries_;
	int * Indices_;
	double * Values_;
	int * TransNumNz_;
	int ** TransIndices_;
	double ** TransValues_;
	int * TransMyGlobalEquations_;
	bool OrigMatrixIsCrsMatrix_;

};

#endif /* EPETRA_CRSMATRIXTRANSPOSER_H */
