//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_GAASP_INTERFACE_H
#define Rythmos_GAASP_INTERFACE_H

#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Rythmos_GAASPErrorEstimate.hpp"
#include "GModelBase.h"
#include "InitGaaspOO.h"
#include "GForwardSolve.h"
#include "GAdjointSolve.h"
#include "GErrorEstimate.h"
#include "GMeshRefine.h"
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

namespace Rythmos {

class GAASPInterface 
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<GAASPInterface> {
  public:
    GAASPInterface();
    ~GAASPInterface() {};
    void forwardSolve();
    void adjointSolve();
    Teuchos::RCP<const GAASPErrorEstimate> computeErrorEstimate();
    void refineMesh();
    void setThyraModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<double> > tModel);
    void initialize(); // Set up GAASP objects after we have a model
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    
    // Redefined from Teuchos::Describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();


  private:

    bool isInitialized_;
    Teuchos::RCP<Thyra::ModelEvaluator<double> > tModel_;
    boost::shared_ptr<GAASP::GModelBase> gModel_;
    double sTime_;    // Start time
    double eTime_;    // End time
    double timeStep_; // Time step size
    GAASP::TMethodFlag method_; // Time integration method pair
    GAASP::TErrorFlag qtyOfInterest_; // Quantity of Interest
    std::vector<double> initialCondition_; // Model DE initial condition
    int dim_;  // Size of solution vector
    GAASP::TRefinementFlag meshRefineMethod_; // Mesh refinement method

    Teuchos::RCP<GAASP::GForwardSolve> forwardSolver_;
    bool forwardSolveCompleted_;
    double **fSoln_;  // Forward solution
    double *mesh_;    // Forward solution mesh
    double **fDer_;   // Forward solution derivatives
    int nSteps_;      // (size of the Forward solution mesh_) - 1

    Teuchos::RCP<GAASP::GAdjointSolve> adjointSolver_;
    bool adjointSolveCompleted_;
    double **aSoln_;   // Adjoint solution

    Teuchos::RCP<GAASP::GErrorEstimate> errorEstimator_;
    bool errorEstimateCompleted_;
    double errorEstimate_; // Total global error estimate
    double **intError_;  // Global Error Approximation by equation by time-step

    Teuchos::RCP<GAASP::GMeshRefine> meshRefiner_;
    bool refineMeshCompleted_;
    double uTOL_;  // Error tolerance for global error (global error control)
    double *newMesh_;    // Refined forward mesh
    int nNodes_;    // Number of nodes in newMesh_

    Teuchos::RCP<Teuchos::ParameterList> paramList_;
    static const std::string sTime_name_;
    static const double sTime_default_;
    static const std::string eTime_name_;
    static const double eTime_default_;
    static const std::string timeStep_name_;
    static const double timeStep_default_;
    static const std::string method_name_;
    static const std::string method_default_;
    static const std::string qtyOfInterest_name_;
    static const std::string qtyOfInterest_default_;
    static const std::string uTOL_name_;
    static const double uTOL_default_;
    static const std::string meshRefineMethod_name_;
    static const std::string meshRefineMethod_default_;


};

} // namespace Rythmos


#endif // Rythmos_GAASP_INTERFACE_H
