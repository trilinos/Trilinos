
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

#ifndef _EPETRA_MUMPS_H_
#define _EPETRA_MUMPS_H_

class Epetra_Import;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_MultiVector;
#include "Epetra_SerialDenseVector.h"
class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseMatrix;
class Amesos_EpetraInterface;
class Amesos_EpetraRedistributor;

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Amesos_EpetraRedistributor.h"

extern "C" {
#include "dmumps_c.h"
}

//! Amesos_Mumps:  An object-oriented wrapper for Mumps.
/*!  Amesos_Mumps will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Mumps solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

  Mumps execution can be tuned through a variety of parameters.
  Amesos_Mumps.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.
    
*/
class Amesos_Mumps : public Amesos_EpetraRedistributor { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Mumps Constructor.
  /*! Creates an Amesos_Mumps instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Mumps(const Epetra_LinearProblem& LinearProblem, const AMESOS::Parameter::List &ParameterList );

  //! Amesos_Mumps Destructor.
  /*! Completely deletes an Amesos_Mumps object.  
  */
  ~Amesos_Mumps(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      </ul>

      postconditions:<ul>
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      <li>MDS will be modified to reflect the symbolic factorization 
      which has been performed.
      </ul>

    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization (and symbolic
      factorization if necessary) on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      </ul>

      postconditions:<ul>
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      <li>MDS will be modified to reflect the numeric factorization 
      which has been performed.
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

#if 0
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Get a pointer to the ParameterList.
  const AMESOS::Parameter::List *GetParameterList() const { return(ParameterList_); };

  //! Returns true if MUMPS can handle this matrix shape 
  /*! Returns true if the matrix shape is one that MUMPS can
    handle. MUMPS only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};


  //! Returns the Schur complement matrix as an Epetra_CrsMatrix.
  /*! Returns the (dense) SSchur complement matrix as an Epetra_CrsMatrix. This
      matrix is defined on all the processes in the Epetra Communicator. However,
      it has rows on the host process only.
      If \in flag : if \c true, MUMPS will compute the Schur complement matrix,
      with respect to the (global) rows defined in the integer array
      \c SchurComplementRows, of size \c NumSchurComplementRows.
      Those two arrays are defined on the host only.
  */
  int ComputeSchurComplement(bool flag,
			     int NumSchurComplementRows, int * SchurComplementRows);

  Epetra_CrsMatrix * GetCrsSchurComplement();

  Epetra_SerialDenseMatrix * GetDenseSchurComplement();
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Reads the parameter list and updates internal variables. 
  /*!
    ReadParameterList is called by SymbolicFactorization.  Hence, few codes 
    will need to make an explicit call to ReadParameterList.
   */
  int ReadParameterList() ;
  //@}

  //! Set prescaling.
  /*! Use double precision vectors of size N (global dimension of the matrix) as
      scaling for columns and rows. \c ColSca and \c RowSca must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(8) = -1.
  */
  int SetPrecscaling(double * ColSca, double * RowSca )
  {
    ColSca_ = ColSca;
    RowSca_ = RowSca;
    return 0;
  }

  //! Set ordering.
  /*! Use integer vectors of size N (global dimension of the matrix) as
      given ordering. \c PermIn must be defined on the host
      only, and allocated by the user, if the user sets ICNTL(7) = 1.
  */
  int SetOrdering(int * PermIn)
  {
    PermIn_ = PermIn;
    return 0;
  }

  int SetMaxis(int Maxis)
  {
    Maxis_ = Maxis;
    return 0;
  }

  int SetMaxs( int Maxs) 
  {
    Maxs_ = Maxs;
    return 0;
  }

  //! Get the pointer to the RINFO array (defined on all processes).
  double * GetRINFO() 
  {
    return (MDS.rinfo);
  }

  //! Get the pointer to the INFO array (defined on all processes).
  int * GetINFO() 
  {
    return (MDS.info);
  }

  //! Get the pointer to the RINFOG array (defined on host only).
  double * GetRINFOG()
  {
    return (MDS.rinfog);
  }

  //! Get the pointer to the INFOG array (defined on host only).
  int * GetINFOG()
  {
    return (MDS.infog);
  }

  int SetKeepMatrixDistributed(bool flag) 
  {
    KeepMatrixDistributed_ = flag;
    return 0;
  }
  
  int SetICNTL(int pos, int value);
  int SetCNTL(int pos, double value);
  
  int PrintInformation();

private:  
  int ConvertToSerial(); 
  /*
    ConvertToTriplet - Convert matirx to form expected by Mumps: Row, Col, Val
    Preconditions:
      numentries_, NumGloalElements_ and SerialMatrix_ must be set.
    Postconditions:
      Row, Col, Val are resized and populated with values that reflect
      the input matrix A.

  */
  int ConvertToTriplet();     

  /*
    PerformSymbolicFactorization - Call Mumps to perform symbolic factorization
    Preconditions:
      Row, Col and Val must be set

    Postconditions:
      MDS will be modified to reflect the symbolic factorization 
      which has been performed.
      SymbolicFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
      
  int PerformSymbolicFactorization(); 

  /*
    PerformNumericFactorization - Call Mumps to perform numeric factorization
    Preconditions:
      IsLocal must be set to 1 if the input matrix is entirely stored on process 0
      Row, Col and Val must be set
    Postconditions:
      MDS will be modified to reflect the numeric factorization 
      which has been performed.
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

  void SetICNTLandCNTL();
  
 protected:

  bool SymbolicFactorizationOK_;   // True if SymbolicFactorization has been done
  bool NumericFactorizationOK_;    // True if NumericFactorization has been done
  bool IsConvertToTripletOK_;
  bool IsComputeSchurComplementOK_;
  
  DMUMPS_STRUC_C MDS ;             // Mumps data structure 

  //
  //  Row, Col, Val form the triplet representation used by Mumps
  //
  Epetra_IntSerialDenseVector * Row; // MS // store COO format in epetra vectors
  Epetra_IntSerialDenseVector * Col;
  Epetra_SerialDenseVector    * Val;

  int MyPID;               //  Process number (i.e. Comm().MyPID() 
  
  int numentries_;         //  Number of non-zero entries in Problem_->GetOperator()
  int NumGlobalElements_;  //  Number of rows and columns in the Problem_->GetOperator()

  bool  KeepMatrixDistributed_;          // MS // this governs the ICNTL(18) parameter.
                                         // MS // If false, then matrix is redistributed
                                         // MS // to proc 0 before converting it to
                                         // MS // triplet format. Then, MUMPS will take care
                                         // MS // of reditribution. If true, the input
                                         // MS // distributed matrix is passed to MUMPS.
  
  const Epetra_Map * Map_;

  int NumMUMPSNonzeros_;                  // MS // actual number of nonzeros in the matrix
  int ErrorMsgLevel_;                     // MS // output level 
  
  bool UseTranspose_;
  
  const AMESOS::Parameter::List * ParameterList_ ; 

  int icntl_[40];                         // MS // to allow users overwrite default settings
  double cntl_[5];                        // MS // as specified by Amesos
  double * RowSca_, * ColSca_;
  int * PermIn_;
  int Maxis_, Maxs_;  

  int NumSchurComplementRows_;            // MS // Schur complement section
  int * SchurComplementRows_;

  Epetra_CrsMatrix * CrsSchurComplement_;
  Epetra_SerialDenseMatrix * DenseSchurComplement_;

  int verbose_;
  
};  // End of  class Amesos_Mumps

#endif /* _EPETRA_MUMPS_H_ */
