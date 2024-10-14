/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef EPETRAEXT_MULTIPOINTMODELEVALUATOR_H
#define EPETRAEXT_MULTIPOINTMODELEVALUATOR_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include "Epetra_ConfigDefs.h"
#include "EpetraExt_ModelEvaluator.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockMultiVector.h"
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif

/** \brief Epetra-based Model Evaluator subclass for Charon!
 *
 * This class will support a wide number of different types of abstract
 * problem types that will allow NOX, LOCA, Rythmos, Aristos, and MOOCHO to
 * solve different types of problems with Charon.
 *
 * ToDo: Finish documentation!
 */

namespace EpetraExt {
  class MultiPointModelEvaluator
    : public ModelEvaluator
  {
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief . */
  MultiPointModelEvaluator(
    Teuchos::RCP<EpetraExt::ModelEvaluator> underlyingME_,
    const Teuchos::RCP<EpetraExt::MultiComm> &globalComm_,
    const std::vector<Epetra_Vector*> initGuessVec,
    Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Vector> > >  qvec,
    Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Vector> > >  matching_vec = Teuchos::null
    );

  //@}

  ~MultiPointModelEvaluator();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

  private:

   Teuchos::RCP<EpetraExt::ModelEvaluator> underlyingME;

   //! Pointer to the global (full XYZT) communicator.
   Teuchos::RCP<EpetraExt::MultiComm> globalComm;

   //! Array of parameter vectors that define the multi-point problem
   Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Vector> > > q_vec;

   //! Pointer to split (spatial) Jacobian matrix.
   Teuchos::RCP<Epetra_RowMatrix> split_W;

   //! Split (spatial) input vector -- local storage.
   Teuchos::RCP<Epetra_Vector> split_x;

   //! Split (spatial) residual vector -- local storage.
   Teuchos::RCP<Epetra_Vector> split_f;

   //! Split vector of response functions -- local storage.
   Teuchos::RCP<Epetra_Vector> split_g;

   //! Split sensitivity vector -- local storage.
   Teuchos::RCP<Epetra_MultiVector> split_DfDp;

   //! Split sensitivity vector -- local storage.
   Teuchos::RCP<Epetra_MultiVector> split_DgDx;
   Teuchos::RCP<Epetra_MultiVector> split_DgDp;

   EpetraExt::ModelEvaluator::DerivativeMultiVector* derivMV_DfDp;
   EpetraExt::ModelEvaluator::Derivative* deriv_DfDp;
   EpetraExt::ModelEvaluator::DerivativeMultiVector* derivMV_DgDx;
   EpetraExt::ModelEvaluator::Derivative* deriv_DgDx;
   EpetraExt::ModelEvaluator::DerivativeMultiVector* derivMV_DgDp;
   EpetraExt::ModelEvaluator::Derivative* deriv_DgDp;

   //! Pointer to global XYZT Jacobian matrix
   Teuchos::RCP<EpetraExt::BlockCrsMatrix> block_W;

   //! Pointer to global multipoint solution vector -- local storage.
   EpetraExt::BlockVector* block_x;

   //! Pointer to global multipoint residual vector -- local storage.
   EpetraExt::BlockVector* block_f;

   //! Pointer to global multipoint DfDp multi vector -- local storage.
   EpetraExt::BlockMultiVector* block_DfDp;

   //! Pointer to global multipoint DfDp multi vector -- local storage.
   EpetraExt::BlockMultiVector* block_DgDx;

   //! Pointer to initial multipoint solution vector.
   Teuchos::RCP<EpetraExt::BlockVector> solution_init;

   //! Number of g vectors supported by underlyingME, often used as a bool
   int underlyingNg;

   //! Number of time steps computed on each time domain.
   int timeStepsOnTimeDomain;

   //! Total number of time step domains.
   int numTimeDomains;

   //! Time domain on current processor.
   int timeDomain;

   /*!
     \brief Stencil for each row of global XYZT Jacobian matrix.

     Used in creating global XYZT Jacobian matrix for different
     finite difference schemes.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   std::vector< std::vector<int> >* rowStencil_int;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   std::vector< std::vector<long long> >* rowStencil_LL;
#endif

   //! Set of indices into global XYZT Jacobian matrix.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   std::vector<int>* rowIndex_int;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   std::vector<long long>* rowIndex_LL;
#endif
   bool longlong;

   //! Some local data
   EDerivativeMultiVectorOrientation orientation_DgDp;
   int num_dg0dp0;
   int num_g0;
   int num_p0;

   //! Array of vectors that have data for g-matching optimization problem
   Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Vector> > > matching_vec;
   bool matchingProblem;

  };
}
#endif
