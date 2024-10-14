// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/*!
 * \file Amesos_Taucs.h
 *
 * \class Amesos_Taucs
 *
 * \brief Interface to TAUCS
 *
 * \date Last updated on 24-May-05.
 */

#ifndef AMESOS_TAUCS_H
#define AMESOS_TAUCS_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Amesos_Control.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
extern "C" {
  //  #include "taucs.h"
}

/*!
  \brief Amesos_Taucs: An interface to the TAUCS package.

  \author Marzio Sala, SNL 9214

  \date Last updated on June 2005
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Taucs_Pimpl ; 
#endif

class Amesos_Taucs: public Amesos_BaseSolver, 
                    private Amesos_Time, 
                    private Amesos_NoCopiable, 
                    private Amesos_Utils, 
                    private Amesos_Control, 
                    private Amesos_Status { 

public: 

  //@{ \name Constructor methods

  //! Default constructor.
  Amesos_Taucs(const Epetra_LinearProblem& LinearProblem );

  //! Default destructor.
  ~Amesos_Taucs(void);
  
  //@}
  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  const Epetra_LinearProblem* GetProblem() const { return(Problem_); };

  //@}
  //@{ \name Query methods.

  bool MatrixShapeOK() const;

  //!  Amesos_Taucs supports only symmetric matrices, hence transpose is irrelevant, but harmless
  /*!  
   */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  bool UseTranspose() const {return(UseTranspose_);};

  const Epetra_Comm& Comm() const {return(GetProblem()->GetOperator()->Comm());};

  int SetParameters( Teuchos::ParameterList &ParameterList);

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( Amesos_Status::NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( Amesos_Status::NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( Amesos_Status::NumSolve_ ); }

  //! Prints timing information
  void PrintTiming() const;
  
  //! Prints status information
  void PrintStatus() const;
 
  //! Extracts timing information from the current solver and places it in the parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming(TimingParameterList); }
 
private:  

  //@}
  //@{ \name Utility Methods

  //! Returns a reference to the RowMatrixRowMap().
  inline const Epetra_Map& Map() const
  {
    return(Matrix_->RowMatrixRowMap());
  }
  
  //! Returns a reference to the linear system matrix.
  inline const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Returns a reference to the already allocated SerialMap.
  inline Epetra_Map& SerialMap() 
  {
    return(*(SerialMap_.get()));
  }
  
  //! Returns a reference to the SerialMatrix.
  inline Epetra_RowMatrix& SerialMatrix()
  {
    return(*(SerialMatrix_.get()));
  }

  //! Returns a reference to the already SerialMatrix as Crs (if allocated).
  inline Epetra_CrsMatrix& SerialCrsMatrix()
  {
    return(*(SerialCrsMatrix_.get()));
  }

  //! Returns a reference to the already allocated Importer.
  inline Epetra_Import& Importer()
  {
    return(*(Importer_.get()));
  }
  
  //! Constructs a matrix with all rows on processor 0.
  int ConvertToSerial();

  //! Converts the Epetra_RowMatrix into TAUCS format.
  int ConvertToTaucs();
  
  //! Performs the symbolic factorization.
  int PerformSymbolicFactorization();

  //! Performs the numeric factorization.
  int PerformNumericFactorization(); 

  // @}
  
  //! If \c true, the transpose of A is used.
  bool UseTranspose_;

  Teuchos::RCP<Epetra_Map> SerialMap_;
  Teuchos::RCP<Epetra_CrsMatrix> SerialCrsMatrix_;
  Teuchos::RCP<Epetra_RowMatrix> SerialMatrix_;
  Teuchos::RCP<Epetra_Import> Importer_;

  const Epetra_Map* Map_;
  const Epetra_RowMatrix* Matrix_;

  //! Pointer to the linear system problem.
  const Epetra_LinearProblem* Problem_;

  //! Quick accessor pointer to internal timing data.
  int MtxConvTime_, MtxRedistTime_, VecRedistTime_;
  int SymFactTime_, NumFactTime_, SolveTime_;

  //
  //  PrivateTaucsData_ contains pointers to data needed by taucs whose
  //  data structures are defined by taucs.h
  //
  Teuchos::RCP<Amesos_Taucs_Pimpl> PrivateTaucsData_; 


};  // class Amesos_Taucs  

#endif
