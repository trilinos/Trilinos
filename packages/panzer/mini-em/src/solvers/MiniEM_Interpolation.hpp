// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_Interpolation_hpp_
#define _MiniEM_Interpolation_hpp_

#include "Teko_Utilities.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestCallback.hpp"

#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif

#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif
#include "MiniEM_Utils.hpp"


/** \brief Registers a mini_em::InterpolationRequestCallback with reqHandler under the given request name, so that the interpolation operator between the lo_basis and ho_basis DOF spaces can later be requested by name.
  *
  * \param[in] name request name to register the interpolation operator under.
  * \param[in] linObjFactory linear object factory used to build the interpolation operator.
  * \param[in] reqHandler the Teko request handler to register the callback with.
  * \param[in] lo_basis_name name of the low-order DOF basis to interpolate from.
  * \param[in] ho_basis_name name of the high-order DOF basis to interpolate to.
  * \param[in] op the Intrepid2 operator (e.g. values, gradients) to interpolate.
  * \param[in] waitForRequest if true, defer building the operator until first requested; if false, build it eagerly.
  * \param[in] dump if true, write the built operator out to a Matrix Market file.
  * \param[in] worksetSize number of cells to process per workset when building the operator.
  * \param[in] matrixFree if true, build a matrix-free interpolation operator instead of an assembled matrix.
  */
void addInterpolationToRequestHandler(const std::string& name,
                                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                                      const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                                      const std::string& lo_basis_name,
                                      const std::string& ho_basis_name,
                                      Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                                      const bool waitForRequest=true,
                                      const bool dump=false,
                                      const size_t worksetSize=1000,
                                      const bool matrixFree=false);


/** \brief Implements the Teko request handler interface to build (on
  * first request, or eagerly if waitForRequest is false) and return
  * an interpolation operator from a low-order to a high-order DOF
  * basis, for use as a preconditioner building block (e.g. for
  * p-multigrid).
  */
class InterpolationRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
private:

  std::string name_;
  const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;
  const std::string lo_basis_name_;
  const std::string ho_basis_name_;
  Intrepid2::EOperator op_;
  const bool dump_;
  Teko::LinearOp interp_;
  const size_t worksetSize_;
  const bool matrixFree_;

public:

  /** \brief Construct from the request name, linear object factory, low-/high-order basis names, the Intrepid2 operator to interpolate (default: values), and options controlling when/how the operator is built.
    *
    * \param[in] name request name this callback handles.
    * \param[in] linObjFactory linear object factory used to build the interpolation operator.
    * \param[in] lo_basis_name name of the low-order DOF basis to interpolate from.
    * \param[in] ho_basis_name name of the high-order DOF basis to interpolate to.
    * \param[in] op the Intrepid2 operator (e.g. values, gradients) to interpolate.
    * \param[in] waitForRequest if true, defer building the operator until first requested; if false, build it eagerly in the constructor.
    * \param[in] dump if true, write the built operator out to a Matrix Market file.
    * \param[in] worksetSize number of cells to process per workset when building the operator.
    * \param[in] matrixFree if true, build a matrix-free interpolation operator instead of an assembled matrix.
    */
  InterpolationRequestCallback(const std::string& name,
                               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                               const std::string& lo_basis_name,
                               const std::string& ho_basis_name,
                               Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                               const bool waitForRequest=true,
                               const bool dump=false,
                               const size_t worksetSize=1000,
                               const bool matrixFree=false);

  /// \brief Builds the interpolation operator (matrix-free or assembled, per the matrixFree option), optionally writing it out to a Matrix Market file if dump is true.
  void build();

  /// \brief Returns true if the requested operator's name matches the name given at construction.
  bool handlesRequest(const Teko::RequestMesg & rm);

  /// \brief Returns the interpolation operator, building it first if it hasn't been built yet.
  Teko::LinearOp request(const Teko::RequestMesg & rm);

  /// \brief Asserts that the requested operator's name matches the name given at construction.
  void preRequest(const Teko::RequestMesg & rm);
};



#endif
