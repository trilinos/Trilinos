/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_LINEARPROBLEMREDISTOR_H
#define EPETRA_LINEARPROBLEMREDISTOR_H

class Epetra_Map; 
class Epetra_Export;
class Epetra_LinearProblem;
class Epetra_RowMatrixTransposer;

//! Epetra_LinearProblemRedistor: A class for redistributing an Epetra_LinearProblem object.

/*! This class provides capabilities to redistribute an existing Epetra_LinearProblem object 
	  across a parallel distributed memory machine.  All or part of a linear problem object can
		be redistributed.   Reverse distributions, value updates and matrix transposition 
		are also supported.  Specification of the redistribution can be done by 
		<ol> 
		<li> providing a target	Epetra_Map object describing the new distribution, or 
		<li> by specifying the number of processors to use and stating whether or not
		     the problem should be completely replicated on all processors.
		</ol>

*/

class Epetra_LinearProblemRedistor {
    
  public:

    //! @name Constructors/destructors
  //@{ 
  //! Epetra_LinearProblemRedistor constructor using pre-defined layout.
  /*!
    \param Problem (In) An existing Epetra_LinearProblem object.  The Epetra_RowMatrix, the LHS and RHS pointers
		       do not need to be defined before this constructor is called.
		\param RedistMap (In) An Epetra_Map describing the target layout of the redistribution.

    \return Pointer to a Epetra_LinearProblemRedistor object.

  */ 
  Epetra_LinearProblemRedistor(Epetra_LinearProblem * OrigProblem, const Epetra_Map & RedistMap);

  //! Epetra_LinearProblemRedistor constructor specifying number of processor and replication bool.
  /*!
    \param Problem (In) An existing Epetra_LinearProblem object.  The Epetra_RowMatrix, the LHS and RHS pointers
		       do not need to be defined before this constructor is called.
		\param NumProc (In) Number of processors to use when redistributing the problem.  Must be between 1 and the
		       number of processor on the parallel machine.
		\param Replicate (In) A bool that indicates if the linear problem should be fully replicated on all processors.
		       If true, then a complete copy of the linear problem will be made on each processor.  If false, then
					 the problem will be roughly evenly spread across the total number of processors.

    \return Pointer to a Epetra_LinearProblemRedistor object.
  */ 
  Epetra_LinearProblemRedistor(Epetra_LinearProblem * OrigProblem, int NumProc, bool Replicate);

  //! Epetra_LinearProblemRedistor copy constructor.
  
  Epetra_LinearProblemRedistor(const Epetra_LinearProblemRedistor& Source);
  
  //! Epetra_LinearProblemRedistor destructor.
  
  virtual ~Epetra_LinearProblemRedistor();
  //@}
  
  //! @name Forward transformation methods
  //@{ 
  
  //! Generate a new Epetra_LinearProblem as a redistribution of the one passed into the constructor.
  /*! Constructs a new Epetra_LinearProblem that is a copy of the one passed in to the constructor.
		  The new problem will have redistributed copies of the RowMatrix, LHS and RHS from the original
			problem.  If any of these three objects are 0 pointers, then the corresponding pointer will be
			zero in the redistributed object.  

			The redistributed matrix will constructed as an Epetra_CrsMatrix.  The LHS and RHS will be Epetra_MultiVector
			objects.

			Two bools can be set when calling this method.  The first,
			ConstructTranspose, will cause the Redistribute method to construct the transpose of the original
			row matrix.  The second, MakeDataContiguous, forces the memory layout of the output matrix, RHS and LHS
			to be stored so that it is compatible with Fortran.  In particular, the Epetra_CrsMatrix is stored so 
			that value from row to row are contiguous, as are the indices.  This is compatible with the Harwell-Boeing
			compressed row and compressed column format.  The RHS and LHS are created so that there is a constant stride between
			the columns of the multivector.  This is compatible with Fortran 2D array storage.

    \param ConstructTranspose (In) Causes the output matrix to be transposed.  This feature can be used
		       to support solvers that need the matrix to be stored in column format.  This option
					 has no impact on the LHS and RHS of the output problem.
		\param MakeDataContiguous (In) Causes the output matrix, LHS and RHS to be stored in a form compatible with
		       Fortran-style solvers.  The output matrix will be compatible with the Harwell-Boeing compressed
					 column format.  The RHS and LHS will be stored such that the last value in column j of the 
					 multivector is stored next to the first value in column j+1.
		\param RedistProblem (Out) The redistributed Epetra_LinearProblem.  The RowMatrix, LHS and RHS that are generated
		       as part of this problem will be destroyed when the Epetra_LinearProblemRedistor object is destroyed.

		\return Integer error code, 0 if no errors, positive value if one or more of the
		        LHS or RHS pointers were 0.  Negative if some other fatal error occured.
					 
  */
  int CreateRedistProblem(const bool ConstructTranspose, const bool MakeDataContiguous, 
													Epetra_LinearProblem *& RedistProblem);

	
  //! Update the values of an already-redistributed problem.
  /*! Updates the values of an already-redistributed problem.  This method allows updating 
		  the redistributed problem without allocating new storage.  All three objects in the RedistProblem will be
			updated, namely the Matrix, LHS and RHS.  If the LHS or RHS are 0 pointers, they will be ignored.

    \param ProblemWithNewValues (In) The values from ProblemWithNewValues will be copied into the RedistProblem.  The
		       ProblemWithNewValues object must be identical in structure to the Epetra_LinearProblem object used to create
					 this instance of Epetra_LinearProblemRedistor.

		\return Integer error code, 0 if no errors, positive value if one or more of the input 
		        LHS or RHS pointers were 0.  Negative if some other fatal error occured.
					 
  */
  int UpdateRedistProblemValues(Epetra_LinearProblem * ProblemWithNewValues);

  //! Update the values of an already-redistributed RHS.
  /*! Updates the values of an already-redistributed RHS.  This method allows updating 
		  the redistributed RHS without allocating new storage.  This method updates only the RHS, and no
			other part of the RedistLinearProblem.

    \param RHSWithNewValues (In) The values from RHSWithNewValues will be copied into the RHS of the RedistProblem.  The
		       RHSWithNewValues object must be identical in structure to the Epetra_MultiVector object used to create
					 this instance of Epetra_LinearProblemRedistor.

		\return Integer error code, 0 if no errors.
					 
  */
  int UpdateRedistRHS(Epetra_MultiVector * RHSWithNewValues);
  //@}
  
  //! @name Reverse transformation methods
  //@{ 
  //! Update LHS of original Linear Problem object.
  /*! Copies the values from the LHS of the RedistProblem Object into the LHS passed in to the method.  If the
		  RedistProblem is replicated, the LHS will be computed as an average from all processor.  
    \param LHS (Out) On exit, the values in LHS will contain the values from the current RedistLinearProblem LHS.
		       If the GIDs of the RedistMap are not one-to-one, e.g., if the map is replicated, the output values for
					 each GID will be an average of all values at that GID.  If the RedistProblem is being solved redundantly
					 in any fashion, this approach to computing the values of LHS should produce a valid answer no matter
					 how the RedistProblem LHS is distributed.
			
    \return Error code, returns 0 if no error.
  */
	int UpdateOriginalLHS(Epetra_MultiVector * LHS);
  //@}
  
  //! @name Attribute accessor methods
  //@{ 
  //! Returns const reference to the Epetra_Map that describes the layout of the RedistLinearProblem.
  /*! The RedistMap object can be used to construct other Epetra_DistObject objects whose maps are compatible with
		  the redistributed linear problem map.
			\warning This method must not be called until \e after CreateRedistProblem() is called.
	*/
	const Epetra_Map & RedistMap() const {return(*RedistMap_);};
  
  //! Returns const reference to the Epetra_Export object used to redistribute the original linear problem.
  /*! The RedistExporter object can be used to redistribute other Epetra_DistObject objects whose maps are compatible with
		  the original linear problem map, or a reverse distribution for objects compatible with the RedistMap().
			\warning This method must not be called until \e after CreateRedistProblem() is called.
  */
  const Epetra_Export & RedistExporter()  const {return(*RedistExporter_);};
  //@}
  
  //! @name Utility methods
  //@{ 

	//! Extract the redistributed problem data in a form usable for other codes that require Harwell-Boeing format.
	/*! This method extract data from the linear problem for use with other packages, such as SuperLU, that require
		  the matrix, rhs and lhs in Harwell-Boeing format.  Note that the arrays returned by this method are owned by the
			Epetra_LinearProblemRedistor class, and they will be deleted when the owning class is destroyed.
	*/
	int ExtractHbData(int & M, int & N, int & nz, int * & ptr, int * & ind, 
										double * & val, int & Nrhs, double * & rhs, int & ldrhs, 
																								double * & lhs, int & ldlhs) const;
  //@}
  
  //! @name I/O methods
  //@{ 
  
  //! Print method
  virtual void Print(ostream & os) const;
  //@}
  
 private: 
	int GenerateRedistMap();

	Epetra_LinearProblem * OrigProblem_;
	int NumProc_;
	Epetra_LinearProblem * RedistProblem_;
	Epetra_Map * RedistMap_;
	Epetra_RowMatrixTransposer * Transposer_;
	Epetra_Export * RedistExporter_;

	bool Replicate_;
	bool ConstructTranspose_;
	bool MakeDataContiguous_;
	bool MapGenerated_;
	bool RedistProblemCreated_;

	mutable int * ptr_;
		

};

#endif /* EPETRA_LINEARPROBLEMREDISTOR_H */
