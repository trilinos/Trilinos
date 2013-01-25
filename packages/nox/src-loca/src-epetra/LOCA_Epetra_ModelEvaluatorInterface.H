// $Id$ 
// $Source$ 

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

#ifndef LOCA_EPETRA_MODEL_EVALUATOR_INTERFACE_H
#define LOCA_EPETRA_MODEL_EVALUATOR_INTERFACE_H

#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_Observer.H"
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "LOCA_DerivUtils.H"
#include "LOCA_Parameter_Vector.H"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "LOCA_Epetra_Group.H"

namespace LOCA {
  namespace Epetra {

    /*!
     * \brief Wrapper for an EpetraExt::ModelEvaluator
     *
     * If an application interfaces their code to solvers using the 
     * EpetraExt::ModelEvaluator, this class provides a wrapper so that the 
     * model evaluator can be used instead of having the user write concrete 
     * versions of the LOCA::Epetra::Interface objects.
     */
    class ModelEvaluatorInterface : 
      public LOCA::Epetra::Interface::TimeDependent,
      public NOX::Epetra::ModelEvaluatorInterface,
      public LOCA::DerivUtils {
    
    public:
    
      //! Constructor
      ModelEvaluatorInterface(
		     const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     const Teuchos::RefCountPtr<EpetraExt::ModelEvaluator>& m,
		     double perturb = 1.0e-6);
    
      //! Destructor
      virtual ~ModelEvaluatorInterface();

      //! Return LOCA parameter vector
      const ParameterVector&  getLOCAParameterVector() const;
    
      virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			    const FillType fillFlag);
    
      virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);
      
      virtual bool computePreconditioner(const Epetra_Vector& x, 
					 Epetra_Operator& M,
					 Teuchos::ParameterList* precParams = 0);

      virtual void setParameters(const ParameterVector& p);

      virtual bool computeShiftedMatrix(double alpha, double beta, 
					const Epetra_Vector& x,
					Epetra_Operator& A);

      virtual void setXdot(const Epetra_Vector& xdot, const double time);

      virtual void
      printSolution  (const Epetra_Vector  &x_, double conParam);

      void
      setObserver( const Teuchos::RCP<NOX::Epetra::Observer> & observer_ );

      // LOCA::DerivUtils interface

      //! Copy constructor
      ModelEvaluatorInterface(const ModelEvaluatorInterface&);

      //! Clone
      virtual Teuchos::RCP<DerivUtils> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      /*! 
       * \brief Compute derivative of f with respect to parameter, identified 
       * by param_id
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeDfDp(LOCA::MultiContinuation::AbstractGroup& grp, 
		  const std::vector<int>& param_ids,
		  NOX::Abstract::MultiVector& result,
		  bool isValidF) const;

      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.
       */
      void
      postProcessContinuationStep(
                           LOCA::Abstract::Iterator::StepStatus stepStatus,
                           LOCA::Epetra::Group& group);
    
    protected:

      //! Parameter vector as an Epetra_Vector
      Epetra_Vector param_vec;

      //! LOCA parameter vector
      LOCA::ParameterVector loca_param_vec;

      //! Epetra vector storing x_dot, filled with zeros
      mutable Epetra_Vector *x_dot;

      //! Stored values of alpha and beta parameters from previous call
      double alpha_prev, beta_prev;

      //! Observer for \c printSolution().
      Teuchos::RCP<NOX::Epetra::Observer> observer;
    
    };

  } // namespace Epetra

} // namespace LOCA

#endif
