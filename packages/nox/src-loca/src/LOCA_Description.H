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

/*! 
\page loca_overview LOCA Overview

\section loca_intro Introduction
  
%LOCA is generic continuation and bifurcation analysis package that is designed for large-scale applications. The algorithms are designed with minimal additional interface requirements over that needed for a Newton method to reach an equilibrium solution. %LOCA is built upon the <a href="./index.html">NOX</a> nonlinear solver package.  The algorithms in %NOX are generic and written to the NOX::Abstract::Group and NOX::Abstract::Vector, which provide abstract interfaces to the linear algebra, data structures, and nonlinear equations to be solved. %LOCA uses the %NOX interface and extends it via additional abstract groups that provide the interface needed for continuation and bifurcation tracking, such as setting parameters and computing derivatives with respect to parameters. 

%LOCA provides several generic groups that take the %NOX group representing the equilibrium equations and implement the extended sets of nonlinear equations representing various forms of continuation and bifurcations (such as the additional equation for arc-length continuation). These extended groups also include generic algorithms for computing the Newton step for the extended system based on the Newton step for the equilibrium equations (e.g. Sherman-Morrison-Woodbury formula for arclength continuation). They are then sent to %NOX for nonlinear solution.  

Finally, %LOCA provides a stepper class that repeatedly calls the %NOX nonlinear solver to compute points along a continuation curve. The design allows for continuation of bifurcations so two-parameter bifurcation sets can be generated. The stepper class relies on several support classes that compute predictors, step sizes, etc.  Each of these are discussed in more detail n the \ref loca_class_overview.

Unlike %NOX which can provide a range of nonlinear solvers using a single abstract interface to the nonlinear equations and linear algebra, %LOCA provides several different levels of functionality, each requiring additional information from the underlying problem.  Therefore, the interface to %LOCA is split among several abstract classes each encapsulating a different level of functionality.  To interface to %LOCA, the user need only provide implementations of those abstract classes for the functionality the user is interested in.  %LOCA provides two complete interfaces:

 - LOCA::LAPACK::Group
 - LOCA::Epetra::Group

both of which implement the interface required for all levels of functionality provided by %LOCA.

\section loca_user LOCA User Information

See \ref loca_user_info.

\section loca_classes Class Overview

See \ref loca_class_overview.

\section loca_params Parameter Specifications

See \ref loca_parameters.

\section loca_tutorials Tutorials

See \ref loca_continuation_tutorial and \ref loca_tp_continuation_tutorial.
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_class_overview LOCA Class Overview

\section continuation Continuation Classes

The lowest level of functionality provided by %LOCA is that for continuation.  The LOCA::MultiContinuation::AbstractGroup provides the interface required for single and multiple parameter continuation and is derived from the NOX::Abstract::Group.  The additional methods provided primarily consist of setting/retrieving parameter values and computing the derivative \f$\partial F/\partial p\f$ for a given parameter \f$p\f$.  For application codes that cannot provide this derivative directly, the LOCA::MultiContinuation::FiniteDifferenceGroup provides a generic implementation using finite differencing.

Multiple continuation algorithms have been implemented, including:

- Natural Continuation (LOCA::MultiContinuation::NaturalGroup)
- Arc-length Continuation (LOCA::MultiContinuation::ArcLengthGroup)

The continuation groups implement the extended set of nonlinear equations to be solved at each continuation step.  

\section bifurcation Bifurcation Classes

%LOCA can locate and track several codimension 1 bifurcations:

- Turning point (a.k.a. fold) (LOCA::TurningPoint::MooreSpence::ExtendedGroup or LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup)
- Pitchfork (LOCA::Pitchfork::MooreSpence::ExtendedGroup)
- Hopf (LOCA::Hopf::MooreSpence::ExtendedGroup).

For each bifurcation, the ExtendedGroup for that bifurcation implements a set of nonlinear equations representing that bifurcation.  Each of these groups provide a complete implementation of the LOCA::MultiContinuation::AbstractGroup and therefore can be used in continuation to track the bifurcation in a second parameter.  The interfaces for these bifurcation groups are encapsulated in the following abstract groups:

- LOCA::TurningPoint::MooreSpence::AbstractGroup (turning points using the Moore-Spence formulation)
- LOCA::TurningPoint::MinimallyAugmented::AbstractGroup (turning points using the minimally augmented formulation)
- LOCA::Pitchfork::MooreSpence::AbstractGroup (pitchfork bifurcations)
- LOCA::Hopf::MooreSpence::AbstractGroup (Hopf bifurcations)

See the relevant class documentation for more details on the additional interfaces provided by these abstract classes.

\section stepping Parameter Stepping Classes

The final component of %LOCA consists of several parameter stepping classes.  These consist of the LOCA::Stepper which drives a %NOX nonlinear solver at each parameter step, and several support classes that compute predictors (classes belonging to the LOCA::MultiPredictor namespace) and encapsulate various step size control strategies (classes in the LOCA::StepSize namespace).  These classes use a strategy/factory design concept, allowing users to provide their own predictors, step-size control algorithms, etc...  Please look at the \ref new_loca_framework page for more information on these classes.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_parameters LOCA Parameter Reference Page

As in %NOX, all parameters are passed to %LOCA through parameter lists defined by the Teuchos::ParameterList class.  All parameters are passed to the LOCA::Stepper class during construction.  Below is a listing of all of the parameters that can be passed to the stepper with links to relevant classes that use those parameters.

\section stepper_params Stepper Parameters

The parameter list passed to the Stepper has two sublists, "LOCA" and "NOX".  The "LOCA" sublist contains parameters and sublists for the Stepper and various %LOCA support classes, while the "NOX" sublist contains parameters and sublists specific to %NOX solvers.  The %NOX sublist is passed directly to an instance of the NOX::Solver::Manager class in the Stepper and therefore has the same structure as any %NOX solver parameter list hierarchy.  Please see the <a href="./parameters.html">%NOX Parameter Reference</a> page for more details on the structure of this list.

<ul> 
<li> <B>"NOX"</B> sublist -- See <a href="./parameters.html">%NOX Parameter Reference</a> for more details
<li> <B>"LOCA"</B> sublist
     <ul>
     <li> <B>"Stepper"</B> sublist -- used by LOCA::Stepper to set parameters relevant for continuation run
          <ul>
          <li> <B>"Continuation Method"</B> -- [string] (default: "Arc Length") Type of continuation to use.  Choices are:
	       <ul>
               <li> "Natural" -- Natural continuation (LOCA::MultiContinuation::NaturalGroup)
               <li> "Arc Length" -- Pseudo arc-length continuation using bordering (LOCA::MultiContinuation::ArcLengthGroup)
                </ul>
	  <li> <B>"Continuation Parameter"</B> -- [string] (Must be supplied) Name of continuation parameter
	  <li> <B>"Initial Value"</B> -- [double] (Must be supplied) Initial value of continuation parameter
	  <li> <B>"Max Value"</B> -- [double] (Must be supplied) Maximum value of continuation parameter
	  <li> <B>"Min Value"</B> -- [double] (Must be supplied) Minimum value of continuation parameter
	  <li> <B>"Max Steps"</B> -- [int] (default: 100) Maximum number of continuation steps (including failed steps)
	  <li> <B>"Max Nonlinear Iterations"</B> -- [int] (default: 15) Maximum number of nonlinear iterations per continuation step
	  <li> <B>"Skip Parameter Derivative"</B> -- [bool] (default: true)  For natural continuation, skip the df/dp computation which is usually unnecessary.
          <li> <B>"Enable Arc Length Scaling"</B> -- [bool] (default: true) Enable arc-length scaling to equilibrate solution and parameter components to arc-length equations (see LOCA::MultiContinuation::ArcLengthGroup)
	  <li> <B>"Goal Arc Length Parameter Contribution"</B> -- [double] (default: 0.5) Goal for parameter contribution to arc-length equation
	  <li> <B>"Max Arc Length Parameter Contribution"</B> -- [double] (default: 0.8) Max for parameter contribution to arc-length equation, triggering rescaling
	  <li> <B>"Initial Scale Factor"</B> -- [double] (default: 1.0) Initial scale factor for parameter term of arc-length equation
	  <li> <B>"Min Scale Factor"</B> -- [double] (default: 1.0e-3) Minimum scale factor for scaling parameter term of arc-length equation
	  <li> <B>"Enable Tangent Factor Step Size Scaling"</B> -- [bool] (default: false) Enable step size scaling by cosine between two consective tangent vectors \f$v_0\f$ and \f$v_1\f$ to continuation curve:  \f$|v_0\cdot v_1|^{\alpha}\f$ where \f$\alpha\f$ is the tangent factor exponent.
	  <li> <B>"Min Tangent Factor"</B> -- [double] (default: 0.1) Minumum cosine between two consecutive tangent vectors, below which the continuation step is failed
	  <li> <B>"Tangent Factor Exponent"</B> -- [double] (default: 1.0) Exponent on the cosine between two consecutive tangent vectors, which then modifies the step size
          <li> <B>"Bordered Solver Method"</B> -- [string] (default:  "Bordering") Method for solving bordered system of equations.  Used by LOCA::BorderedSolver::Factory.  Available choices are
	        <ul>
		<li> "Bordering" -- Bordering (LOCA::BorderedSolver::Bordering)
                <li> "Nested" -- Nested for method recursive bordered solves (LOCA::BorderedSolver::Nested).  This method is appropriate when the system is itself bordered, for example, when using the minimally augmented turning point formulation or when the system has constraints.  This method combines the nested bordered rows and columns and forms one large bordered system.  The solver used for the resulting bordered system is specified through the "Nested Bordered Solver" sublist.
                      <ul>
		      <li> <B>"Nested Bordered Solver"</B> sublist -- Specifies the solver for the resulting bordered system in a nested bordered system.  Any valid bordered solver parameter may be used here.
		      </ul>
                </ul>
	  <li> <B>"Compute Eigenvalues"</B> -- [bool] (default: false) Flag for requesting eigenvalue calculation after each continuation step
          <li> <B>"Eigensolver"</B> sublist -- used by LOCA::Eigensolver::Factory to determine eigensolver strategies
          <ul> 
          <li> <B>"Method"</B> -- [string] (default "Default") Choices are
                <ul>
                <li> "Default" -- No eigenvalues are computed (LOCA::Eigensolver::DefaultStrategy)
                <li> "Anasazi" -- Eigenvalues/vectors are computed using %Anasazi (LOCA::Eigensolver::AnasaziStrategy).  %Anasazi parameters are
                    <ul>
                    <li> <B>"Operator"</B> -- [string] (default: "Jacobian inverse") Operator to compute eigenvalues of.   Used by LOCA::AnasaziOperator::Factory.  Valid choices are:
                         <ul>
                         <li> "Jacobian Inverse" -- Eigenvalues of \f$J^{-1}\f$ (LOCA::AnasaziOperator::JacobianInverse).
                         <li> "Shift-Invert" -- Eigenvalues of \f$(J-\sigma M)^{-1}M\f$ (LOCA::AnasaziOperator::ShiftInvert).
                             <ul>
                             <li> <B>"Shift"</B> -- [double] (default 0.0) \f$\sigma\f$ as defined above
                             </ul>
                         <li> "Cayley" -- Eigenvalues of \f$(J-\sigma M)^{-1}(J - \mu M)\f$ (LOCA::AnasaziOperator::Cayley).
                             <ul>
                             <li> <B>"Cayley Pole"</B> -- [double] (default 0.0) \f$\sigma\f$ as defined above
                             <li> <B>"Cayley Zero"</B> -- [double] (default 0.0) \f$\mu\f$ as defined above
                             </ul>
                         </ul>
	            <li> <B>"Block Size"</B> -- [int] (default: 1) Block size
	            <li> <B>"Num Blocks"</B> -- [int] (default: 30) Maximum number of blocks (equals the maximum length of the Arnoldi factorization)
	            <li> <B>"Num Eigenvalues"</B> -- [int] (default: 4) Number of requested eigenvalues
	            <li> <B>"Convergence Tolerance"</B> -- [double] (default: 1.0e-7) Tolerance for the converged eigenvalues
	            <li> <B>"Linear Solve Tolerance"</B> -- [double] (default: uses the same tolerance from the nonlinear solver) Tolerance for the linear solves used by the eigensolver.
	            <li> <B>"Step Size"</B> -- [int] (default: 1) Checks convergence every so many steps
	            <li> <B>"Maximum Restarts"</B> -- [int] (default: 1) Number of restarts allowed
                    <li> <B>"Symmetric"</B> -- [bool] (default:  false) Is the operator symmetric
                    <li> <B>"Normalize Eigenvectors with Mass Matrix"</B> -- [bool] (default:  false) Option to normalize vMv=1 instead of vv=1.
	            <li> <B>"Verbosity"</B> -- [Anasazi::MsgType] (default: Anasazi::Errors+Anasazi::Warnings+Anasazi::FinalSummary) Verbosity level
	            <li> <B>"Sorting Order"</B> -- [string] (default: "LM") Sorting order of printed eigenvalues.  Used by LOCA::EigenvalueSort::Factory.  Valid choices are:
		         <ul>
			 <li> "LM" -- Largest magnitude (LOCA::EigenvalueSort::LargestMagnitude)
			 <li> "LR" -- Largest real component (LOCA::EigenvalueSort::LargestReal)
			 <li> "LI" -- Largest imaginary component (LOCA::EigenvalueSort::LargestImaginary)
			 <li> "SM" -- Smallest magnitude (LOCA::EigenvalueSort::SmallestMagnitude)
			 <li> "SR" -- Smallest real component (LOCA::EigenvalueSort::SmallestReal)
			 <li> "SI" -- Smallest imaginary component (LOCA::EigenvalueSort::SmallestImaginary)
			 <li> "CA" -- Largest real part of inverse Cayley transformation (LOCA::EigenvalueSort::LargestRealInverseCayley)
			 </ul>
                    </ul>
	        </ul>
          </ul>
	  </ul>
     <li> <B>"Bifurcation"</B> sublist -- used by LOCA::Bifurcation::Factory to determine what type of bifurcation calculation, if any, to use
          <ul>
	  <li> <B>"Bifurcation Parameter"</B> -- [string] (Must be supplied if <B>"Type"</B> is not <B>"None"</B>) Name of bifurcation parameter
          <li> <B>"Type"</B> -- [string] (default: "None") Bifurcation method to use.  Choices are
               <ul>
               <li> "None" -- No bifurcation and do equilibrium continuation instead
               <li> "Turning Point" -- Turning point bifurcations
	            <ul> 
		    <li> <B>"Formulation"</B> -- Name of the bifurcation formulation. Valid choices are
		         <ul>
                         <li> "Moore-Spence" -- [Default] Moore-Spence formulation (LOCA::TurningPoint::MooreSpence::ExtendedGroup)
			      <ul>
			      <li> <B>"Length Normalization Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Pointer to NOX::Abstract::Vector holding length normalization vector.
			      <li> <B>"Initial Null Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Initial guess for null vector.
			      <li> <B>"Perturb Initial Solution"</B> -- [bool] (default: false) Flag indicating whether to apply an initial perturbation to the initial guess for the solution vector before starting bifurcation algorithm
			      <li> <B>"Relative Perturbation Size"</B> -- [double] (default: 1.0e-3) Size of relative perturbation of initial guess for solution vector.
			      <li> <B>"Solver Method"</B> -- [string] (default:  "Salinger Bordering") Solver method for solving extended turning point equations.  Valid choices are
			          <ul>
				  <li> "Salinger Bordering" -- Original bordering algorithm that uses (nearly) singular Jacobian solves (LOCA::TurningPoint::MooreSpence::SalingerBordering).
				  <li> "Phipps Bordering" -- Bordering algorithm that uses nonsingular bordered matrix solves (LOCA::TurningPoint::MooreSpence::PhippsBordering).
				       <ul>
				       <li> <B>"Bordered Solver Method"</B> -- Method for solving bordered systems of equations.  See description in the "Stepper" sublist above.
				       </ul>
				  </ul>
			      </ul>
                         <li> "Minimally Augmented" -- Minmally augmented formulation (LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup)
			      <ul>
			      <li> <B>"Bordered Solver Method"</B> -- Method for solving bordered systems of equations.  See description in the "Stepper" sublist above.
			      <li> <B>"Symmetric Jacobian"</B>  -- [bool] (default: false) Flag indicating whether Jacobian matrix is symmetric.
			      <li> <B>"Constraint Method"</B> -- [string] (default: "Default") Type of constraint method to use.  Valid choices are
			          <ul>
				  <li> "Default" -- Default method (LOCA::TurningPoint::MinimallyAugmented::Constraint)
				  <li> "Modified" -- A modified method that computes updates to the null vectors every nonlinear interation, instead of directly solving for them (LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint)
				  </ul>
			      <li> <B>"Initial Null Vector Compuation"</B> -- [string] (default: "User Provided") - Method to compute initial \f$a\f$ and \f$b\f$ vectors.  Valid choices are:
			           <ul>
				   <li> "User Provided" - Initial vectors are provided in the parameter list, in which case the following parameters are relevant:
				        <ul>
					<li> <B>"Initial A Vector"</B>  -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Vector storing initial value for \f$a\f$ vector
					<li> <B>"Initial B Vector"</B>  -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied for nonsymmetric Jacobians) Vector storing initial value for \f$b\f$ vector
					</ul>
				   <li> "Solve df/dp" -- Compute \f$a = J^{-T}df/dp\f$ and \f$b = J^{-1} df/dp\f$ where \f$p\f$ is the bifurcation parameter.
				   </ul>
			      <li> <B>"Update Null Vectors Every Continuation Step"</B>  -- [bool] (default: true) Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors continuation step
			      <li> <B>"Update Null Vectors Every Nonlinear Iteration"</B>  -- [bool] (default: false) Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors every nonlinear iteration
			      </ul>
			 </ul>
		    </ul>
               <li> <B>"Pitchfork"</B> -- Pitchfork bifurcations
	            <ul>
	            <li> <B>"Formulation"</B> -- Name of the bifurcation formulation. Valid choices are
		         <ul>
		         <li> "Moore-Spence" [Default] -- Moore-Spence formulation (LOCA::Pitchfork::MooreSpence::ExtendedGroup)
			      <ul>
			      <li> <B>"Length Normalization Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Pointer to NOX::Abstract::Vector holding length normalization vector.
			      <li> <B>"Initial Null Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Initial guess for null vector.
			      <li> <B>"Asymmetric Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Pointer to NOX::Abstract::Vector holding asymmetric vector \f$\psi\f$ used in pitchfork algorithm.
			      <li> <B>"Perturb Initial Solution"</B> -- [bool] (default: false) Flag indicating whether to apply an initial perturbation to the initial guess for the solution vector before starting bifurcation algorithm
			      <li> <B>"Relative Perturbation Size"</B> -- [double] (default: 1.0e-3) Size of relative perturbation of initial guess for solution vector.
			      <li> <B>"Solver Method"</B> -- [string] (default:  "Salinger Bordering") Solver method for solving extended turning point equations.  Valid choices are
			          <ul>
				  <li> "Salinger Bordering" -- Original bordering algorithm that uses (nearly) singular Jacobian solves (LOCA::Pitchfork::MooreSpence::SalingerBordering).
				  <li> "Phipps Bordering" -- Bordering algorithm that uses nonsingular bordered matrix solves (LOCA::Pitchfork::MooreSpence::PhippsBordering).
				       <ul>
				       <li> <B>"Bordered Solver Method"</B> -- Method for solving bordered systems of equations.  See description in the "Stepper" sublist above.
				       </ul>
				  </ul>
			      </ul>
			 <li> "Minimally Augmented" -- Minmally augmented formulation (LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup)
			      <ul>
			      <li> <B>"Bordered Solver Method"</B> -- Method for solving bordered systems of equations.  See description in the "Stepper" sublist above.
			      <li> <B>"Symmetric Jacobian"</B>  -- [bool] (default: false) Flag indicating whether Jacobian matrix is symmetric.
			      <li> <B>"Asymmetric Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Pointer to NOX::Abstract::Vector holding asymmetric vector \f$\psi\f$ used in pitchfork algorithm.
			      <li> <B>"Initial Null Vector Compuation"</B> -- [string] (default: "User Provided") - Method to compute initial \f$a\f$ and \f$b\f$ vectors.  Valid choices are:
			           <ul>
				   <li> "User Provided" - Initial vectors are provided in the parameter list, in which case the following parameters are relevant:
				        <ul>
					<li> <B>"Initial A Vector"</B>  -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Vector storing initial value for \f$a\f$ vector
					<li> <B>"Initial B Vector"</B>  -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied for nonsymmetric Jacobians) Vector storing initial value for \f$b\f$ vector
					</ul>
				   <li> "Solve df/dp" -- Compute \f$a = J^{-T}df/dp\f$ and \f$b = J^{-1} df/dp\f$ where \f$p\f$ is the bifurcation parameter.
				   </ul>
			      <li> <B>"Update Null Vectors Every Continuation Step"</B>  -- [bool] (default: true) Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors continuation step
			      <li> <B>"Update Null Vectors Every Nonlinear Iteration"</B>  -- [bool] (default: false) Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors every nonlinear iteration
			      </ul>
                         </ul>
		    </ul>
               <li> <B>"Hopf"</B> -- Hopf bifurcations
	            <ul>
	            <li> <B>"Formulation"</B> -- Name of the bifurcation formulation. Valid choices are
		         <ul>
			 <li> "Moore-Spence" -- [Default] -- Moore-Spence formulation (LOCA::Hopf::MooreSpence::ExtendedGroup)
			      <ul>
			      <li> <B>"Length Normalization Vector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Pointer to NOX::Abstract::Vector holding length normalization vector.
			      <li> <B>"Initial Real Eigenvector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Initial guess for real component of generalized eigenvector in Hopf algorithm
			      <li> <B>"Initial Imaginary  Eigenvector"</B> -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Initial guess for imaginary component of generalized eigenvector in Hopf algorithm
			      <li> <B>"Initial Frequency"</B> -- [double] (Must be supplied) Initial guess for Hopf frequency \f$\omega\f$ in Hopf algorithm
			      <li> <B>"Perturb Initial Solution"</B> -- [bool] (default: false) Flag indicating whether to apply an initial perturbation to the initial guess for the solution vector before starting bifurcation algorithm
			      <li> <B>"Relative Perturbation Size"</B> -- [double] (default: 1.0e-3) Size of relative perturbation of initial guess for solution vector.
			      <li> <B>"Solver Method"</B> -- [string] (default:  "Salinger Bordering") Solver method for solving extended turning point equations.  Valid choices are
			          <ul>
				  <li> "Salinger Bordering" -- Original bordering algorithm that uses (nearly) singular matrix solves (LOCA::Hopf::MooreSpence::SalingerBordering).
				  </ul>
			      </ul>
			 </ul>
		    </ul>
	       <li> <B>"User-Defined"</B> User defined bifurcation.  Set 
                   <B>"User-Defined Name"</B> to be the parameter list name containing user-defined strategy, which must be of type Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>.
               </ul>
          </ul>
     <li> <B>"Predictor"</B> sublist -- used by LOCA::MultiPredictor::Factory to determine what type of predictor to use for each continuation step
          <ul>
	  <li> <B>"Method"</B> -- [string] (default: "Secant") Predictor method to use for computing the initial guess for each continuation step.  Choices are 
                 <ul>
                 <li> "Constant" -- Use previous solution as initial guess for next step (LOCA::MultiPredictor::Constant)
                 <li> "Secant" -- Use secant vector to continuation curve to compute initial guess (LOCA::MultiPredictor::Secant)
                 <li> "Tangent" -- Use tangent vector to continuation curve to compute initial guess (LOCA::MultiPredictor::Tangent)
                 <li> "Random" -- Use a random perturbation of previous solution to compute initial guess (LOCA::MultiPredictor::Random)
                 </ul> 
	  <li> <B>"Epsilon"</B> -- [double] (default: 1.0e-3) Relative size of perturbation for random predictor
	  <li> <B>"First Step Predictor"</B> sublist -- used by the secant predictor to determine which predictor to use for the first continuation step -- Replicates <B>"Predictor"</B> sublist structure
	  <li> <B>"Last Step Predictor"</B> sublist -- Predictor to use for last step of arc-length continuation to hit target (max or min) value exactly (usually "Constant" or "Random") -- Replicates <B>"Predictor"</B> sublist structure
	  </ul>
     <li> <B>"Step Size"</B> sublist -- used by LOCA::StepSize::Factory to determine step size constrol strategies
          <ul>
	  <li> <B>"Method"</B> -- [string] (default: "Adaptive") Step size control strategy to use.  Choices are
                <ul>
                <li> "Constant" -- Use a constant step size in general, reducing the step size after a failure and increasing step size back up to original value after subsequent successes (see LOCA::StepSize::Constant)
                <li> "Adaptive" -- Use an adaptive step size control strategy that adjusts step size according to the number of Newton iterations per step (see LOCA::StepSize::Adaptive)
                </ul>
	  <li> <B>"Initial Step Size"</B> -- [double] (default: 1.0) Initial parameter step size
	  <li> <B>"Min Step Size"</B> -- [double] (default: 1.0e-12) Minimum parameter step size
	  <li> <B>"Max Step Size"</B> -- [double] (default: 1.0e+12) Maximum parameter step size
	  <li> <B>"Failed Step Reduction Factor"</B> -- [double] (default: 0.5) Factor by which step size is reduced after a failed step
	  <li> <B>"Successful Step Increase Factor"</B> - [double] (default: 1.26) Factor by which step size is increased after a successful step when the step size is smaller than the initial step size (Constant step size method only). 
	  <li> <B>"Aggressiveness"</B> -- [double] (default: 0.5) Aggressiveness factor in adaptive step size adjustment
	  </ul>
      <li> <B>"Constraints"</B> sublist -- Used to provide additional constraint equations.
          <ul>
	  <li> <B>"Constraint Object"</B> -- [Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>] (Must be supplied) Pointer to object implementig constraint interface LOCA::MultiContinuation::ConstraintInterface.
	  <li> <B>"Constraint Parameter Names"</B> -- [Teuchos::RCP< std::vector<std::string> >] (Must be supplied) Pointer to vector of string names of parameters that become additional unkowns for solving constrained system.
	  <li> <B>"Bordered Solver Method"</B> -- Method for solving bordered systems of equations arising from constraints.  See description in the "Stepper" sublist above.
	  </ul>
      </ul>
</ul>

\section user_defined_list User-Defined Strategies Through Parameter Lists

Most strategy choices above have another available option that is not listed above for providing user-defined strategies.  To use a user-defined strategy, one must implement the corresponding strategy abstract interface (e.g., LOCA::MultiPredictor::AbstractStrategy) and instantiate an object of the user-defined type.  Then set the corrsponding "Method" (or similar) parameter of the relevant sublist to "User-Defined", set the parameter "User-Defined Name" to the parameter list name of the strategy (e.g., "My Predictor"), and set the parameter with that name to be a ref-count pointer to that object with the type of the abstract strategy interface (e.g., Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>).  

\section user_defined_factory User-Defined Strategies Through Factories

Another option is available for providing user-defined strategies by implementing a factory to instantiate these strategies.  This option is most appropriate when multiple user-defined strategies are available, and the user wishes to control which strategy is instantiated through parameter lists.   This is often used when providing strategies that are specific to the underlying %NOX group, for example the %Epetra and %LAPACK factories described below.  In this case one imlements the LOCA::Abstract::Factory by defining the relevant strategy creation methods and then passes a ref-count pointer to this factory in the LOCA::createGlobalData() call.  Then the user can select these strategies just as if they were built in to %LOCA.

\subsection epetra_factory Epetra Factory

The LOCA::Epetra::Factory provides strategies specific to the NOX::Epetra::Group.  It provides the following strategies:
<ul>
<li> Methods for solving bordered systems of equations (LOCA::BorderedSolver::Factory)
    <ul>
    <li> <B>"Householder"</B> -- Use the QR (a.k.a. Householder) method (LOCA::BorderedSolver::EpetraHouseholder).
       <ul>
       <li> <B>"Include UV In Preconditioner"</B> -- [bool] (default: false) Flag indicating whether to use the \f$U\f$ and \f$V\f$ terms in the preconditioner for \f$P\f$
       <li> <B>"Use P For Preconditioner"</B> -- [bool] (default: false)  Flag indicating whether to use the representation of \f$P\f$ as a LOCA::Epetra::LowRankUpdateRowMatrix for computing the preconditioner.  This is valid only for preconditioners that accept an Epetra_RowMatrix interface.
       <li> <B>"Transpose Solver Method"</B> -- [string] (default: "Transpose Preconditioner") Method for preconditioning the transpose linear system (LOCA::Epetra::TransposeLinearSystem::Factory).  Available choices are:
       <ul>
	  <li> "Transpose Preconditioner" -- Use the transpose of the preconditioner for the original system.
	  <li> "Left Preconditioning" -- Use the transpose of the preconditioner, and apply using left preconditioning.
	  <li> "Explicit Transpose" -- Form the transpose of the matrix and compute the preconditioner.  This method is available only if Trilinos is configured with %EpetraExt support (--enable-epetraext).
       </ul>
       <li> <B>"Preconditioner Method"</B> -- [string] (default: "Jacobian") - 
           Method for preconditioning the \f$P\f$ operator.  Choices are:
           <ul>
           <li> "Jacobian" (default) -- Use the preconditioner for \f$J\f$
           <li> "SMW" -- Use the Sherman-Morrison-Woodbury formula for the
                inverse of \f$P\f$, replacing the inverse of \f$J\f$ with
                the preconditioner for \f$J\f$.
           </ul>
       <li> <B>"Scale Augmented Rows"</B> -- [bool] (default: true) -
           Scale augmented rows to unit 2-norm before computing QR
           factorization.
      </ul>
    <li> <B>"Augmented"</B> -- Form an augmented linear problem implemented in a pseudo-matrix free mode (LOCA::BorderedSolver::EpetraAugmented).
    </ul>
</ul>

\subsection lapack_factory LAPACK Factory

The LOCA::LAPACK::Factory provides strategies specific to the NOX::LAPACK::Group.  It provides the following strategies:
<ul>
<li> Methods for solving bordered systems of equations (LOCA::BorderedSolver::Factory)
    <ul>
    <li> <B>"LAPACK Direct Solve"</B> -- Form the bordered system as a %LAPACK matrix and solve via %LAPACK.
    </ul>
<li> Eigensolver methods (LOCA::Eigensolver::Factory)
    <ul>
    <li> <B>"DGGEV"</B> -- Use LAPACK's generalized eigenvalue solver DGGEV
       <ul>
       <li> <B>"Num Eigenvalues"</B> -- [int] (default: 4) Number of eigenvalues to print
       <li> <B>"Sorting Order"</B> -- [string] (default: "LM") Sorting order of printed eigenvalues.  Used by LOCA::EigenvalueSort::Factory.  See Eigensolver sublist above for a list of available options.
       </ul>
    </ul>
</ul>
       
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_user_info LOCA User Information

\section loca_user_info_overview Overview

Since %LOCA is built upon %NOX, using %LOCA is quite similar to using %NOX and
an understanding of how to use %NOX is a prerequisite for using %LOCA.  See
\ref nox_user_information for a description of how to use %NOX.  Here we will
discuss the additional steps required to use %LOCA.

\section loca_user_info_step_1 Step 1: Download, Configure, Compile, and Install LOCA

%LOCA is part of %NOX, so by downloading %NOX you have %LOCA.  Note that %LOCA 
is dependent on some Trilinos packages (namely %Teuchos), and greatly benefits
from others (in particular, %Anasazi).  Therefore %LOCA can only be configured
and built inside of Trilinos.  To enable %LOCA 
during configuration and compilation, supply the --enable-loca configuration 
option.  %LOCA %LAPACK and %Epetra support and examples are automatically 
enabled if they are enabled in %NOX (see \ref nox_configuration_options). In
addition, Trilinos must be configured with --enable-teuchos, and if %LAPACK
support is enabled, --enable-teuchos-complex.  %LOCA is capable of using the 
%Anasazi eigensolver to compute eigenvalues/eigenvectors at each continuation
step, and this can be enabled by adding the --enable-anasazi and 
--with-loca-anasazi configure flags. See \ref nox_configuration_options for
other available configure options.  Typing "make" and "make install" in the top-level %NOX directory will then
compile and install %LOCA as well as %NOX.

\section loca_user_info_step_2 Step 2: Create concrete implementations of the LOCA AbstractGroup classes

See \ref loca_class_overview for a description of the %LOCA AbstractGroup class
hierarchy used to supply interfaces for various parameter continuation and
bifurcation tracking algorithms.  These AbstractGroups extended the interface
supplied by the NOX::Abstract::Group to provided the functionality needed for
continuation and bifurcation tracking, and to use %LOCA, concrete 
implementations of these interfaces must be provided.  

\subsection loca_user_info_step_2_writing_implementation Writing your own implementation

As in %NOX, we recommend using the LOCA::LAPACK::Group as a guide for writing
your own implementation.  The simplest approach is to derive your group
from the LOCA::Abstract::Group and provide implementations of the appropriate
AbstractGroup methods for the level of functionality you desire.  At a 
minimum for parameter continuation, you must provide implementations of the
parameter setting/retrieval methods.  Note there is no additional %LOCA 
abstract vector interface to implement to use %LOCA.

\subsection loca_user_info_step_2_loca_implementations Implementations provided with LOCA

%LOCA provides two ready-made implementations:

<ul>
<li>
The LOCA::LAPACK implementation extends the NOX::LAPACK implementation to
%LOCA.  %LAPACK support and examples are automatically enabled if they are
enabled in %NOX.  The LOCA::LAPACK implementation supports all levels of
functionality provided by %LOCA and can be used in any continuation or
bifurcation calculation.

<li>
The LOCA::Epetra implementation extends the NOX::Epetra implementation to
%LOCA.  %Epetra support and examples are automatically enabled if they are
enabled in %NOX.  All levels of functionality provided by LOCA are supported, 
including recently added Hopf and generalized eigenvalue support.

<li>
Currently, there is no %LOCA Petsc support provided, but this will be
added as soon as there is demand.
</ul>

\section loca_user_info_step_3 Step 3: Call LOCA from your code

See \ref loca_continuation_tutorial for a tutorial describing how to use 
%LOCA to perform a simple continuation. The
example problems in the examples-lapack and examples-epetra
directories can serve as templates as well.

\section loca_user_info_step_4 Step 4: Link your code to LOCA

Linking with %LOCA is very similar to linking with %NOX.  You must additionally
link against -lloca and one of the concrete implementations, either 
-llocalapack (if using the LOCA::LAPACK), -llocaepetra (if using 
LOCA::Epetra),or a custom implementation.

\section loca_user_info_problems If you have problems...

If you have problems with %LOCA, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_continuation_tutorial LOCA Continuation Tutorial

\section loca_continuation_tutorial_overview Overview

Here we provide a brief tutorial for using %LOCA to perform a simple
continuation using the Chan problem (see ChanProblemInterface).  The 
code fragements discussed below can be found in ChanContinuation.C in
the Chan subdirectory of the %LOCA %LAPACK examples directory.  

The ChanProblemInterface implements a 1-D finite-difference discretization
of the Chan problem:
\f[
      \frac{d^2 T}{d x^2} + \alpha s \left(1 + \frac{x + 0.5 x^2}{1 + 0.01 x^2}\right) = 0
  \f]
subject to the boundary conditions \f$T(0) = T(1) = \beta\f$.  The parameters
are \f$\alpha\f$, \f$\beta\f$, \f$s\f$, and \f$n\f$, the size of the 
discretization.  The scaling factor \f$s\f$ is used to test continuation
at different scales of \f$\alpha\f$.  Below we will track the solution
\f$T(x,\alpha,\beta)\f$ as \f$\alpha\f$ is varied from \f$0\f$ to \f$5\f$.

\section loca_continuation_tutorial_code ChanContinuation.C Line by Line

\code
#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ChanProblemInterface.H"
\endcode

To use %LOCA in your code, you must always include the header LOCA.H.  Since
this is a %LAPACK problem, we also include LOCA_LAPACK.H.  Finally, we include
the header for the Chan problem, ChanProblemInterface.H.

\code
int main()
{
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;

  alpha = alpha / scale;

  try {
\endcode

Next we set up the basic problem parameters.  \em maxNewtonIters is the
maximum number of nonlinear solver iterations we wish to take at each
continuation step.  %NOX and %LOCA do throw exceptions (but only when
serious errors occur), so all %NOX and %LOCA calculations should be placed
in a try block.

\code
    // Create output file to save solutions
    ofstream outFile("ChanContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;
\endcode

Here we set up a file for saving the solutions computed at each continuation
step.  The printSolution method of the ChanProblemInterface is set up in
such a way that if a file is provided, the current continuation parameter
and solution vector are appended to the file, in addition to printing a portion
of the solution vector to the screen.  The format of this file is a series
of rows, with each row containing \f$n+1\f$ numbers, the first is the 
continuation parameter with the remaining \f$n\f$ consisting of each 
component of the solution vector.

\code
    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");// Default
    //stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Parameter", "alpha");  // Must set
    stepperList.set("Initial Value", alpha);             // Must set
    stepperList.set("Max Value", 5.0/scale);             // Must set
    stepperList.set("Min Value", 0.0/scale);             // Must set
    stepperList.set("Max Steps", 50);                    // Should set
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters); // Should set
    stepperList.set("Compute Eigenvalues",false);        // Default

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");                 // Default

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");               // Default
    //predictorList.set("Method", "Constant");
    //predictorList.set("Method", "Tangent");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");             // Default
    stepSizeList.set("Initial Step Size", 0.1/scale);   // Should set
    stepSizeList.set("Min Step Size", 1.0e-3/scale);    // Should set
    stepSizeList.set("Max Step Size", 10.0/scale);      // Should set
\endcode

Next we set up the %LOCA parameters.  We are setting up the problem to 
perform arc-length continuation in the parameter "alpha" from 0 to 5 with 
a maximum of 50 continuation steps and \em maxNewtonIters nonlinear iterations
per step.  Since we are doing an
equilibrium continuation, we set the bifurcation method to "None".  We use a secant predictor
and adaptive step size control with an initial step size of 0.1, maximum of 
10.0 and minimum of 0.001. 

\code
    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
		      NOX::Utils::Details +
		      NOX::Utils::OuterIteration + 
		      NOX::Utils::InnerIteration + 
		      NOX::Utils::Warning + 
		      NOX::Utils::StepperIteration +
		      NOX::Utils::StepperDetails +
		      NOX::Utils::StepperParameters);  // Should set
\endcode

Next we set up the %NOX parameters.  We use a simple full-step Newton method 
for the nonlinear solve at each continuation step (the default), and specify
what output we want printed to the screen.

\code

    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

\endcode

We now create the LAPACK factory and the global data object.  The factory is optional, and allows the use of LAPACK-specific strategies (e.g., the LAPACK eigensolver DGGEV).  However the global data object is manditory.  If no factory is provided, the second argument to createGlobalData is empty.

\code
    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);
\endcode

Next we instantiate the ChanProblemInterface and create a parameter vector
to store the values of the problem parameters.  Note that it is not 
necessary to put every problem parameter into the parameter vector, only those
that serve as possible continuation or bifurcation parameters need to be
supplied.

\code
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp = 
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));
    
    grp->setParams(p);
\endcode

Next we instantiate the %LAPACK group with the Chan problem and then set
the parameter vector in the group.  From this point on, the LOCA::Stepper,
via the LOCA::LAPACK::Group, will take care of setting parameters in
the problem interface.

\code
    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> normF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Generic> comboOR = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      normF, 
					      maxIters));
\endcode

Next we create appropriate status tests for the problem.  For convergence
at each continuation step, we require the extended (solution and parameter
components) residual norm be smaller than 1.0e-8 and the number of nonlinear
iterations be smaller than \em maxNewtonIters.

\code
    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    // Check for convergence
    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }
\endcode

Finally we instantiate the stepper, run the continuation, and check the
returned status.

\code
    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::LAPACK::Group> finalGroup = 
      Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution = 
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup->getX());

    // Print final solution
    globalData->locaUtils->out()
	            << std::endl << "Final solution is " << std::endl;
    finalGroup->printSolution(finalSolution,
		    	      finalGroup->getParam("alpha"));

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    outFile.close();

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
\endcode

Lastly we copy the final solution out of the stepper, print out the final
parameter list, close the output file, destroy the global data object, and 
catch any thrown exceptions.

After running the example and plotting the maximum of the temperature versus the continuation parameter \f$\alpha\f$ at each step, we  obtain the following continuation curve with two turning points:

\image html chan.png

\image latex chan.eps

Turning point bifurcations occur near \f$\alpha=4\f$ and \f$\alpha=3\f$.  For a tutorial on locating these turning point bifurcations and tracking them in the second parameter \f$\beta\f$, see the \ref loca_tp_continuation_tutorial.

The above plot was generated via MATLAB using the output file ChanContinuation.dat specified above.  For those interested, the MATLAB commands used to generate this plot are shown below.
\code

% open output file
fid = fopen('ChanContinuation.dat');

% read dimension of discretization
n = fscanf(fid, '%d', 1);
  
alpha = []; % array of continuation parameter values at each step
x = [];     % array of solution components at each step
  
% read values from output file
while ~feof(fid)
  
  % read alpha
  alpha = [alpha fscanf(fid, '%g', 1)];
  
  % read x
  x = [x fscanf(fid, '%g', n)];
  
end

% close output file
fclose(fid);

% compute maximum of each temperature profile
maxT = max(x);

plot(alpha,maxT,'bo-');
xlabel('\alpha');
ylabel('T_{max}   ','Rotation',0);
title('Arc-length Continuation:  \beta = 0');

\endcode
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_tp_continuation_tutorial LOCA Turning Point Continuation Tutorial

\section loca_tp_continuation_tutorial_overview Overview

In the \ref loca_continuation_tutorial, two turning point bifurcations were discoved in the Chan problem as the parameter \f$\alpha\f$ was varied from \f$0\f$ to \f$5\f$.  Here we provide a brief tutorial on locating one of these turning points and tracking it in the second parameter \f$\beta\f$.  There are several methods in %LOCA for computing turning points.  The example here uses the minimally augmented formulation using the %LAPACK direct solver for solving the augmented linear equations.  See LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup for a description of this algorithm.  Also see LOCA::TurningPoint::MooreSpence::ExtendedGroup for the Moore-Spence turning point formulation.  The code fragements discussed below can be found in ChanTPContinuation.C in the Chan subdirectory of the %LOCA %LAPACK examples directory.

\section loca_tp_continuation_tutorial_code ChanTPContinuation.C

Much of the setup for the turning point continuation problem is the same as for the equilibrium continuation discussed in \ref loca_continuation_tutorial.  Therefore we will only discuss the differences between the two setups here.

\code
#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ChanProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 4.0;
    double beta = 0.0;
    double scale = 1.0;
    int maxNewtonIters = 10;

    // Create output file to save solutions
    ofstream outFile("ChanTPContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);
\endcode

By examining the plot in \ref loca_continuation_tutorial, a turning point bifurcation occurs near \f$\alpha=4\f$ for \f$\beta = 0\f$.  We use these values as an initial set of parameter values near the bifurcation.  We also set up an output file to store the continuation parameter, solution vector, bifurcation parameter, right null vector, measure of singularity of the Jacobian (sigma) and left null vector at each continuation step.  The format is the same as in ChanContinuation.dat, consisting of a series of rows each containing \f$n+1\f$ numbers.  Three rows are written for each continuation step, the first containing the continuation parameter and the \f$n\f$ components of the solution vector, the second containing the bifurcation parameter and the \f$n\f$ components of the right null vector, and the third containing sigma and the \f$n\f$ components of the left null vector.

\code
    // Create output file to save solutions
    ofstream outFile("ChanTPContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;

    // Create initial guess for the null vector of jacobian
    Teuchos::RCP<NOX::Abstract::Vector> nullVec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    nullVec->init(1.0);               // initial value 1.0

    // Create initial values for a and b for minimally augmented method
    Teuchos::RCP<NOX::Abstract::Vector> a_vec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    a_vec->init(1.0);

    Teuchos::RCP<NOX::Abstract::Vector> b_vec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    b_vec->init(1.0);
\endcode

The only additional set up required for turning point tracking in this problem is to compute an initial guess for the null vector of the Jacobian (for the Moore-Spence formulation) or the initial values for the \f$a\f$ and \f$b\f$ vectors in the minimally augmented formulation.  Here we use a vector of all one's in all three cases.  

\code
    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");  // Default
    stepperList.set("Continuation Parameter", "beta");     // Must set
    stepperList.set("Initial Value", beta);                // Must set
    stepperList.set("Max Value", 1.0);                     // Must set
    stepperList.set("Min Value", 0.0);                     // Must set
    stepperList.set("Max Steps", 20);                      // Should set
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters); // Should set

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Turning Point");          // For turning point
    bifurcationList.set("Bifurcation Parameter", "alpha"); // Must set

    // For Moore-Spence formulation w/bordering
    //bifurcationList.set("Formulation", "Moore-Spence");          // Default
    //bifurcationList.set("Solver Method", "Salinger Bordering");  // Default
    //bifurcationList.set("Solver Method", "Phipps Bordering");
    //bifurcationList.set("Bordered Solver Method",  
    //                    "LAPACK Direct Solve");   // For Phipps Bordering
    //bifurcationList.set("Length Normalization Vector", nullVec); // Must set
    //bifurcationList.set("Initial Null Vector", nullVec);         // Must set

    // For minimally augmented formulation
    bifurcationList.set("Formulation", "Minimally Augmented");
    bifurcationList.set("Initial A Vector", a_vec);                // Must set
    bifurcationList.set("Initial B Vector", b_vec);                // Must set

    // For minimally augmented method, should set these for good performance
    // Direct solve of bordered equations
    bifurcationList.set("Bordered Solver Method",  "LAPACK Direct Solve");
    // Combine arc-length and turning point bordered rows & columns
    stepperList.set("Bordered Solver Method", "Nested");
    Teuchos::ParameterList& nestedList = 
      stepperList.sublist("Nested Bordered Solver");
    // Direct solve of combined bordered system
    nestedList.set("Bordered Solver Method", "LAPACK Direct Solve");
\endcode

We now set \f$\beta\f$ to be the continuation parameter and \f$\alpha\f$ to be the bifurcation parameter.  We will vary \f$\beta\f$ from \f$0\f$ to \f$1\f$, computing a value of \f$\alpha\f$ for each corresponding value of \f$\beta\f$.  The initial value of \f$\alpha\f$ is set internally by accessing the component "alpha" in the parameter vector \em p set below.  In the bifurcation sublist, we indicate that we would like to do turning point tracking using the minimally augmented formulation, and pass ref-count pointers to the initial guess for the a and b vectors.  Note that these must be casted to NOX::Abstract::Vector ref-count pointers.

Both arc-length continuation and the minimally augmented turning point method add one additional equation to be solved resulting in two nested bordered systems, each adding an additional row and column to the Jacobian matrix.  We tell %LOCA to combine these rows and columns into one bordered system with two augmented rows and columns (by setting the "Bordered Solver Method" of the "Stepper" sublist to "Nested") and instruct %LOCA to use the LAPACK-specific linear solver for solving this 2-bordered system.  Note that we must set the "LAPACK Direct Solve" choice twice , once in the "Bifurcation" sublist and once in the "Nested Bordered Solver" sublist.  The first specifies the solver for the first and last continuation steps which use natural continuation instead of arc-length continuation.  The last specifies the solver for the rest of the continuation steps.

\code
     // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");                         // Default

    // Should use for Moore-Spence w/Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& firstStepPredictor 
    //  = predictorList.sublist("First Step Predictor");
    //firstStepPredictor.set("Method", "Random");
    //firstStepPredictor.set("Epsilon", 1.0e-3);

    // Should use for Moore-Spence w/Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& lastStepPredictor 
    //  = predictorList.sublist("Last Step Predictor");
    //lastStepPredictor.set("Method", "Random");
    //lastStepPredictor.set("Epsilon", 1.0e-3);
\endcode

We again use a secant predictor to compute an initial guess at each continuation step.  Because we are using the secant predictor, a different predictor must be used for the first step.  The default is the "Constant" predictor which is fine for the minimally augmented formulation.  However for the Moore-Spence formulation with Salinger bordering, we should use the "Random" predictor because the algorithm is singular.  If the problem Jacobian does not depend on the continuation parameter, we can obtain highly ill-conditioned linear solves when not using the random predictor.  A random predictor can be chosen for the last step for the same reason.

\code
    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");                      // Default
    stepSizeList.set("Initial Step Size", 0.1);                  // Should set
    stepSizeList.set("Min Step Size", 1.0e-3);                   // Should set
    stepSizeList.set("Max Step Size", 1.0);                      // Should set

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
		      NOX::Utils::OuterIteration + 
		      NOX::Utils::OuterIterationStatusTest + 
		      NOX::Utils::InnerIteration +
		      NOX::Utils::Details + 
		      NOX::Utils::Warning + 
		      NOX::Utils::StepperIteration +
		      NOX::Utils::StepperDetails +
		      NOX::Utils::StepperParameters);            // Should set

    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grp = 
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));
    
    grp->setParams(p);

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10, 
					      NOX::StatusTest::NormF::Scaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      statusTestA, statusTestB));

    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    outFile.close();

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
\endcode

The rest of the driver setup is very similar to ChanContinuation.C  We slightly tighten the convergence tolerance to demonstrate the well conditioning of the minimally augmented method.  In this case we must use the LAPACK factory since we are using the LAPACK-specific solver method "LAPACK Direct Solve".

After running the example and reading the data file ChanTPContinuation.dat, we can plot the continuation parameter \f$\beta\f$ versus the bifurcation parameter \f$\alpha\f$ to get a locus of turning point bifurcations in the \f$(\beta,\alpha)\f$ parameter space:

\image html chan_tp.png

\image latex chan_tp.eps

There are two branches of the bifurcation curve which come together to form a cusp.  Starting at \f$\beta=0\f$ on one branch, traversing the cusp, and moving to \f$\beta=0\f$ on the second branch connects the two turning points shown in the \ref loca_continuation_tutorial.

As in the continuation tutorial, the MATLAB code used to generate this plot is shown below.

\code
% open output file
fid = fopen('ChanTPContinuation.dat');

% read dimension of discretization
n = fscanf(fid, '%d', 1);
  
beta = [];     % array of continuation parameter values at each step
alpha = [];    % array of bifurcation parameter values at each step
sigma = [];    % array of singular estimates at each step
x = [];        % array of solution components at each step
y = [];        % array of right null vector components at each step
z = [];        % array of left null vector components at each step

while ~feof(fid)
  
  % read beta
  beta = [beta fscanf(fid, '%g', 1)];
  
  % read x
  x = [x fscanf(fid, '%g', n)];
  
  % read alpha
  alpha = [alpha fscanf(fid, '%g', 1)];
  
  % read y
  y = [y fscanf(fid, '%g', n)];

  % read sigma
  sigma = [sigma fscanf(fid, '%g', 1)];

  % read z
  z = [z fscanf(fid, '%g', n)];
  
end

% close output file
fclose(fid);

% compute maximum of each temperature profile
maxT = max(x);

figure;
plot(beta,alpha,'bo-');
xlabel('\beta');
ylabel('\alpha    ','Rotation',0);
title('Locus of Turning Points');
\endcode
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page new_loca_framework New LOCA Framework

Most of %LOCA has been rewritten to more fully support multi-parameter continuation, constraint tracking, higher-order predictors and advanced bifurcation algorithms.  The new framework completely replaces the old framework, with the old one now removed.  However the interfaces for both the new and old frameworks are still very similar, so the conversion between them is fairly simple.  There are three major components to this refactoring:
<ul>
<li> Use of Teuchos::RCP objects instead of raw pointers or references for better C++ memory management.
<li> Use of strategy interfaces and factories to choose between alternative numerical algorithms.
<li> Use of multi-vectors throughout %LOCA to consolidate single and multiple parameter continuation code and to more fully support block algorithms.
</ul>

In the new framework, Teuchos::RCP's are used whenever objects are passed into an object method and the supplied arguments have a persisting relationship with that object beyond the lifetime of the function call.  This greatly simplifies memory management since the owning object does not need to worry about deleting any of the supplied arguments.  For the most part this is transparent to the user since the only place where users must use ref-count pointers is when passing the group, parameter list, and status test objects to the stepper constructor (or when implementing a strategy interface as discussed below).

In the old framework, whenever a numerical continuation or bifurcation algorithm needed the result of a computation that was problem or data structure dependent, a group method was created for that computation which the user could implement in their own concrete group.  The problem with this however is it ties the implementation of this method to the group and makes it difficult for the group to support multiple ways of implementing the computation.  In the new framework, a strategy interface approach is taken instead where an abstract interface class is written the encapsulates the interface for the computation (strategy interface) and strategy objects are written that implement that interface for different strategies.  Then a factory is written that instantiates strategy objects based on a name provided in a parameter list (see for example LOCA::Eigensolver::AbstractStrategy and LOCA::Eigensolver::Factory or LOCA::BorderedSystem::AbstractStrategy and LOCA::BorderedSystem::Factory).  This allows the same group to use multiple strategies for a given computation and separates the interface to parallel linear algebra data structures in the group from the often many ways data structure dependent computations can be written in the strategy objects.  (In many ways this is similar to the Manager concept in %NOX, except it allows users to easily provide their own strategy implementations by implementing their own factory, see LOCA::Factory and LOCA::Abstract::Factory).

Finally, multi-vectors are used where appropriate in the new framework to consolidate single and multi-parameter continuation code by viewing single parameter continuation as a special case of multi-parameter continuation.  This allows all of the same code to be used for both types of continuation, greatly simplifying code maintenance and development.  It also gives %LOCA a natural and efficient interface to block solvers such as %Anasazi and Belos.
*/
