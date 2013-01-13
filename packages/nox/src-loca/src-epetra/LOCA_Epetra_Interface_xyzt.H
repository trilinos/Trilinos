//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_INTERFACE_XYZT_H
#define LOCA_EPETRA_INTERFACE_XYZT_H

#include "LOCA.H"
#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "LOCA_Epetra_xyztPrec.H"

#ifdef HAVE_NOX_EPETRAEXT

#include <iostream>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DistObject.h"
#include "Epetra_Import.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "EpetraExt_MultiMpiComm.h"
#else
#include "Epetra_SerialComm.h"
#include "EpetraExt_MultiSerialComm.h"
#endif

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"

// Forward declarations
namespace LOCA {
  namespace Epetra {
    namespace Interface {
      class Required;
    }
  }
}
namespace NOX {
  namespace Epetra {
    namespace Interface {
      class Jacobian;
    }
  }
}


namespace LOCA {
  namespace Epetra {
    namespace Interface {
      
      /*!
       * \brief Used by LOCA::Epetra::Group to provide a link to the
       * external code for creating and solving space-time formulations 
       * of the user's problem.
       *
       * Code that takes a standard NOX-LOCA problem interface for 
       * spatially-discretized problems, and creates a NOX-LOCA interface 
       * for space-time problems.  This interface is called by multiple 
       * replicas of the spatially discretized problem for parallelism of 
       * the time domain. Each replica can also own one or more time steps. 
       * This interface assembles a single space-time system including the 
       * solution and residual BlockVectors and the Jacobian BlockCrsMatrix.
       * This interface relies heavily on the EpetraExt_Block* classes. 
       * 
       */
      class xyzt : 
	public LOCA::Epetra::Interface::Required,
	public NOX::Epetra::Interface::Jacobian {
	
      public:
	
	/*! 
	  \brief Constructor

	  The interface, vector, and matrix arguments are all for the spatially
	  decomposed problem. The MultiMpiCOmm communicator gives this
	  class the information needed to construct the space-time system. 
	*/
	xyzt( const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> 
	     &interface,
	     const Epetra_MultiVector &splitMultiVec,
	     const Teuchos::RCP<Epetra_RowMatrix> &splitJac,
	     const Teuchos::RCP<EpetraExt::MultiComm> &globalComm,
	     const Epetra_Vector &initialCondVec,
             double dt,
	     Teuchos::ParameterList *precPrintParams = 0, 
	     Teuchos::ParameterList *precLSParams = 0);
	
	//! Destructor
	virtual ~xyzt();
	
	/*! 
	  \brief NOX Interface function for computing the residual.

	  Computes \em F given the vector \em x. 
	  Returns true if computation was successful.
	*/
	virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			      const FillType fillFlag);
	
	/*! 
	  \brief NOX Interface function for computing the Jacobian matrix. 

	  Computes \em Jac given the vector \em x.
	  Returns true if computation was successful.
	*/
	virtual bool computeJacobian(const Epetra_Vector& x,
				     Epetra_Operator& Jac);
	
	//! LOCA Interface function for setting the continuation parameter
	virtual void setParameters(const LOCA::ParameterVector& param);
	
	//! LOCA Interface function for printing the solution
	virtual void printSolution(const Epetra_Vector& x_, double conParam);
	
	//! Method to alert object that Jacobian fills need to be altered
        //! to compute operator for Floquet theory (monodromy operator)
	void setFloquetFillFlag(bool fff);

	//! Method to start the application of the monodromy matrix,
        //! which involves multiplication by one mass matrix block and
	void beginFloquetOperatorApplication(Epetra_Vector& v);

	//! Method to start the application of the monodromy matrix,
        //! which involves multiplication by one mass matrix block and
	void finishFloquetOperatorApplication(Epetra_Vector& v);
	
	//! Accessor for the global solution BlockVector
	virtual EpetraExt::BlockVector& getSolution();
	
	//! Accessor for the global BlockCrsMatrix
	virtual EpetraExt::BlockCrsMatrix& getJacobian();
	
	//! Accessor for the xyzt Preconditioner
	virtual LOCA::Epetra::xyztPrec& getPreconditioner();
	
	//! LOCA Interface exception handler
	virtual void throwError(const std::string& functionName, 
				const std::string& errorMsg) const;
	
      private:
	
	//! Pointer to interface for usual time-dependent problem
	Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> interface;
	
	//! Pointer to split (spatial) Jacobian matrix.
	Teuchos::RCP<Epetra_RowMatrix> splitJac;
	
	//! Pointer to the global (full XYZT) communicator.
	Teuchos::RCP<EpetraExt::MultiComm> globalComm;
	
	//! Split (spatial) input vector.
	Epetra_Vector splitVec;
	
	//! Split (spatial) residual vector.
	Epetra_Vector splitRes;
	
	//! Split (spatial) input vector from previous step.
	Epetra_Vector splitVecOld; //tmp space
	
	//! Split (spatial) vector with initial condition of transient runs
	const Epetra_Vector initialCondVec;
	
	//! Pointer to global XYZT Jacobian matrix
	EpetraExt::BlockCrsMatrix* jacobian;
	
	//! Pointer to global XYZT solution vector.
	EpetraExt::BlockVector* solution;
	
	//! Pointer to global XYZT solution overlap vector.
	EpetraExt::BlockVector* solutionOverlap;
		
	//! Pointer to global XYZT solution overlap vector importer.
	Epetra_Import* overlapImporter;
	
	//! Number of time steps computed on each time domain.
	int timeStepsOnTimeDomain; 
	
	//! Total number of time step domains.
	int numTimeDomains;
	
	//! Time domain on current processor.
	int timeDomain;
	
	//! Continuation step number (used in printing)
	int conStep;
	
	/*! 
	  \brief Stencil for each row of global XYZT Jacobian matrix.

	  Used in creating global XYZT Jacobian matrix for different 
	  finite difference schemes.
	*/
	std::vector< std::vector<int> >* rowStencil;
	
	//! Set of indices into global XYZT Jacobian matrix.
	std::vector<int>* rowIndex;
	
	/*! 
	  \brief Reference to parameter list for controlling output 
	  of Jacobian matrix right preconditioning steps.
	*/
	Teuchos::ParameterList* precPrintParams;
	
	/*! 
	  \brief Reference to parameter list for controlling solution of
	  linear systems used in Jacobian matrix right preconditioning.
	*/
	Teuchos::ParameterList* precLSParams;
	
	//! Pointer to XYZT preconditioner operator.
	Teuchos::RCP<LOCA::Epetra::xyztPrec> preconditioner;
	
	//! Pointer to split (spatial) Jacobian matrix in CRS matrix form.
	Epetra_CrsMatrix* splitJacCrs;
	
	//! Pointer to split (spatial) mass matrix saved for Floquet operation.
	Epetra_RowMatrix* savedSplitMassForFloquet;
	
	//! Flag indicating whether space matrix is Crs (determined by cast attempt)
	bool  isCrsMatrix;

	//! Flag indicating whether XYZT system is periodic
	bool  isPeriodic;
	
	//! Flag indicating Jacobian fill needs to be altered for Floquet calcs
	bool floquetFillFlag;

	//! Time step size for computing time derivative -- generalize later
        double dt;


	
      };

    } // namespace Interface
  } // namespace Epetra
} // namespace LOCA

#endif // HAVE_NOX_EPETRAEXT
#endif // LOCA_EPETRA_INTERFACE_XYZT_H
