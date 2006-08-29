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
          string        name = ""                 ) ;
  
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
  Teuchos::RefCountPtr<Epetra_Vector> getExactSolution()
  { return exactSolution; }

  // Additional status output
  virtual void outputStatus( ostream & os );


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

  Teuchos::RefCountPtr<Epetra_Import> ColumnToOverlapImporter        ;
  Teuchos::RefCountPtr<Epetra_Vector> oldSolution                    ;
  Teuchos::RefCountPtr<Epetra_Vector> exactSolution                  ;
  Teuchos::RefCountPtr<Epetra_Vector> dTdx                           ;

  int             interface_elem                 ;
  int             local_node                     ;
  int             interface_node                 ;
  int             opposite_node                  ;

  ConvDiff_PDE *  depProbPtr                     ;
};
#endif
