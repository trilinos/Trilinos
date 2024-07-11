// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Finite Element Problem Class
/* Provides function (F) and Jacobian evaluations for the following equation
 * via a 1D linear finite element discretization with Epetra objects.
 *
 * d2u
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
 */

#ifndef _NOX_TEST_CONVDIFF_PDE_EQUATION_H
#define _NOX_TEST_CONVDIFF_PDE_EQUATION_H

#include "GenericEpetraProblem.H" // Base class

// Single PDE Equation Problem Class for HMX Cook-off problem

class ConvDiff_PDE : public GenericEpetraProblem
{

public:

  enum PAPER_SCENARIO {
                        RVALUE_0_00     ,
                        RVALUE_0_30     ,
                        RVALUE_1_00     ,
                        RVALUE_2_50     ,
                        RVALUE_5_67     ,
                        NONE
                      };

  enum INTERFACE_LOCATION {
                                LEFT    ,
                                RIGHT
                          };

  // Constructor
  ConvDiff_PDE(
          Epetra_Comm & Comm                    ,
          double        peclet                  ,
          double        radiation               ,
          double        kappa                   ,
          double        bcWeight                ,
      double        xmin                    ,
      double        xmax                    ,
      double        Tleft                   ,
      double        Tright                  ,
          int           NumGlobalUnknowns = 0   ,
          std::string        name = ""                 ) ;

  // Destructor
  ~ConvDiff_PDE();

  // Method to compute exact interface temperature
  static double computeAnalyticInterfaceTemp(
                  double Radiation      ,
                  double T_left  = 0.0  ,
                  double T_right = 1.0  ,
                  double kappa   = 0.1  ,
                  double peclet  = 9.0      );


  // Set initial condition for solution vector
  void initializeSolution(double val = 0.0);

  // Compute and output analytic slution here
  void initialize();

  // Evaluates the function (F) and/or the Jacobian using the solution
  // values in solnVector.
  bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
                const Epetra_Vector *solnVector,
        Epetra_Vector *rhsVector);

  // Inherited method that allows us to explcitly provide off-block structural info
  // in the form of global composite equation indices
  bool suppliesOffBlockStructure()
  { return expandJacobian; }

  // Inherited method that allows us to explcitly provide off-block structural info
  // in the form of global composite equation indices
  bool suppliesOffBlockStructure(int depID)
  { return expandJacobian; }

  // Method to toggle expanded Jacobian option
  void setExpandJacobian( bool setting )
  { expandJacobian = setting; }

  // Method that allows us to provide off-block structural info
  void getOffBlockIndices( map<int, std::vector<int> > & );

  // Inherited method that allows us to compute fluxes before transferring data
  void prepare_data_for_transfer();

  // Compute the appropriate end heat flux
  void computeHeatFlux( const Epetra_Vector * pVec = NULL );

  // Inherited method that allows us to tidy up transferred data
  void process_transferred_data();

  // Inherited method that allows us to specialize transfers across an interface
  void doTransfer();

  // Accessor
  double getHeatFlux()
  { return myFlux; }

  // Accessor
  double getInterfaceTemp();

  // Accessor
  Teuchos::RCP<Epetra_Vector> getExactSolution()
  { return exactSolution; }

  // Additional status output
  virtual void outputStatus( std::ostream & os );


protected:

  // inserts the global column indices into the Graph
  void generateGraph();

protected:

  double xmin                           ;
  double xmax                           ;
  double Tleft                          ;
  double Tright                         ;
  double dx                             ;
  double T1_exact                       ;

  double peclet                         ;
  double radiation                      ;
  double kappa                          ;

  PAPER_SCENARIO     paperScenario      ;
  INTERFACE_LOCATION myInterface        ;
  double dirScale                       ;
  double myFlux                         ;

  double bcWeight                       ;

  Teuchos::RCP<Epetra_Import> ColumnToOverlapImporter        ;
  Teuchos::RCP<Epetra_Vector> oldSolution                    ;
  Teuchos::RCP<Epetra_Vector> exactSolution                  ;
  Teuchos::RCP<Epetra_Vector> dTdx                           ;

  int             interface_elem        ;
  int             local_node            ;
  int             interface_node        ;
  int             opposite_node         ;

  // Flag to enable off-block structure in compositeJacobian
  bool expandJacobian                   ;

  // Pointer to the other (assumed only one) problem
  ConvDiff_PDE *  depProbPtr            ;
};
#endif
