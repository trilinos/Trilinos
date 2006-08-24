//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_EPETRA_XYZTPREC_H
#define LOCA_EPETRA_XYZTPREC_H

#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Vector.H"
#include "EpetraExt_MultiMpiComm.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

namespace LOCA {
  namespace Epetra {
    
    /*! 
      \brief Preconditioner operator class for solving space-time (XYZT) systems.

      Implements right preconditioning operators for use in global 
      XYZT Jacobian matrix solves. 

      <ul>
        <li> \em Global - applies a right preconditioner to the global 
	      XYZT Jacobian matrix

        <li> \em Sequential - applies single block right
	      preconditioning sequentially in time. This
	      preconditioner is intended as an efficient competitor to
	      the \em Global preconditioner by preconditioning using
	      only the nonzero blocks.

        <li> \em Parallel - simultaneously applies sequential right
              preconditioning across the decoupled time domains. This
              means there is no communication of solutions between
              time doamins.

        <li> \em BlockDiagonal - similar to the \em Parallel
              preconditioner, simultaneously applies sequential right
              preconditioning across the decoupled time domains, but
              only using the diagonal blocks of the Jacobian
              matrix. Note that the BlockDiagonal and Parallel
              preconditioners are equivalent when each time domain
              contains only one time step.

        <li> \em Parareal - two pass right preconditioning applying
              \em Sequential preconditioner over a coarse time grid
              (first time steps of each time domain) and then \em
              Parallel preconditioning across the decoupled time
              domains. This can be thought of as a linearized \em
              parareal strategy for acheiving parallelism in time. The
              benefit over \em Parallel preconditioning alone is that
              an estimate of the solution from the time step on the
              previous domain is computed to help accelerate
              convergence.

        <li> \em BDSDT (block diagonal in space, diagonal in time)

      </ul>
    */
    class xyztPrec : public Epetra_Operator, 
		     public NOX::Epetra::Interface::Required,
		     public NOX::Epetra::Interface::Jacobian,
		     public NOX::Epetra::Interface::Preconditioner
    {
      
    public:
      
      /*!
	\brief Constructor

        Builds a preconditioner operator for a full XYZT Jacobian
	matrix \em jacobian. Right preconditioner applies are
	controlled using the parameters in \em precLSParams.
      */
      xyztPrec(EpetraExt::BlockCrsMatrix& jacobian, 
	       Epetra_CrsMatrix &splitJac,
	       Epetra_CrsMatrix &splitMass,
	       EpetraExt::BlockVector& solution,
	       EpetraExt::BlockVector& solutionOverlap,
	       Epetra_Import &overlapImporter,
	       Teuchos::ParameterList& precPrintParams, 
	       Teuchos::ParameterList& precLSParams, 
	       const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> globalComm_);
      
      //! Destructor
      virtual ~xyztPrec();
      
      //! Set transpose.
      virtual int SetUseTranspose(bool UseTranspose);
      
      //! Apply XYZT preconditioner operator
      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
      //! Apply XYZT preconditioner operator inverse
      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
      
      //! Computing infinity norm
      virtual double NormInf() const;
      
      //! Label
      virtual const char* Label() const;
      
      //! Transpose
      virtual bool UseTranspose() const;
      
      //! Have norm-inf
      virtual bool HasNormInf() const;
      
      /*! 
       * \brief Returns a pointer to the Epetra_Comm communicator associated 
       * with this operator.
       */
      virtual const Epetra_Comm& Comm() const;
      
      /*! 
       * \brief Returns the Epetra_Map object associated with the domain of 
       * this operator.
       */
      virtual const Epetra_Map& OperatorDomainMap() const;
      
      /*! 
       * \brief Returns the Epetra_Map object associated with the range of 
       * this operator.
       */
      virtual const Epetra_Map& OperatorRangeMap() const;
      
      //! Compute residual \f$F\f$.
      virtual bool computeF(const Epetra_Vector&, 
			    Epetra_Vector&, 
			    NOX::Epetra::Interface::Required::FillType);
      
      //! Compute Jacobian \f$J\f$.
      virtual bool computeJacobian(const Epetra_Vector&, Epetra_Operator&);
      
      //! Compute preconditioner \f$M\f$.
      virtual bool computePreconditioner(const Epetra_Vector& x,
					 Epetra_Operator& Prec,
					 Teuchos::ParameterList* p = 0);
      
      
      //! Exception handler for the XYZT preconditioner class.
      virtual void throwError(const string& functionName, 
			      const string& errorMsg) const;
      
      
    private:
      
      //! Reference to global XYZT Jacobian matrix.
      EpetraExt::BlockCrsMatrix& jacobian;
      
      //! Reference to split (spatial) Jacobian matrix.
      Epetra_CrsMatrix &splitJac;
      
      //! Reference to split (spatial) mass matrix.
      Epetra_CrsMatrix &splitMass;
      
      //! Reference to global XYZT solution.
      EpetraExt::BlockVector& solution;
      
      //! Reference to global XYZT solution overlap stored on other processors.
      EpetraExt::BlockVector& solutionOverlap;
      
      //! Reference to Epetra importer for overlapping XYZT solution imports.
      Epetra_Import& overlapImporter;
      
      //! Reference to parameter list for controlling output of solve.
      Teuchos::ParameterList& printParams;
      
      //! Reference to parameter list for controlling linear system solve.
      Teuchos::ParameterList& lsParams;
      
      //! Reference to global communicator for global XYZT operator.
      const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> globalComm;
      
      //! Pointer to global XYZT linear system to be solved.
      std::vector<NOX::Epetra::LinearSystemAztecOO *> linSys;
      
      //! Pointer to jacobian blocks of global XYZT system.
      std::vector<Teuchos::RefCountPtr<Epetra_CrsMatrix> > jacobianBlock;
      
      //! Pointer to mass matrix blocks of global XYZT system.
      std::vector<Teuchos::RefCountPtr<Epetra_CrsMatrix> > massBlock;
      
      //! Pointer to subdiagonal blocks used in the BDSDT perconditioner.
      std::vector<Teuchos::RefCountPtr<Epetra_Vector> > diagBlockSubdiag;
      
      //! Pointer to global XYZT residual vector.
      EpetraExt::BlockVector *residual;
      
      //! Pointer to split (spatial) input vector of current iteration.
      Epetra_Vector *splitVec;
      
      //! Pointer to split (spatial) residual vector of current iteration.
      Epetra_Vector *splitRes;
      
      //! Pointer to split (spatial) input vector from previous iteration.
      Epetra_Vector *splitVecOld;
      
      //! View of split residual vector.
      NOX::Epetra::Vector *splitRes_NEV;
      
      //! View of split input vector.
      NOX::Epetra::Vector *splitVec_NEV;
      
      //! Flag indicating whether XYZT system is periodic.
      bool  isPeriodic;
      
      //! Label for this operator.
      string label;
      
    };
    
  } // namespace Epetra
} // namespace LOCA

#endif // LOCA_EPETRA_XYZTPREC_H
