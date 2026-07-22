// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VALID_PARAMETERS_H
#define ROL_VALID_PARAMETERS_H

#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

/* ROL Parameters */

inline ROL::Ptr<const ROL::ParameterList> getValidROLParameters() {
  
  typedef ROL::ParameterList PL;

  ROL::Ptr<PL> rol = ROL::makePtr<PL>("ROL");

  /* ===== GENERAL INPUT PARAMETERS ============================ */ 
  PL &general = rol->sublist("General");
    general.set("Recompute Objective Function",             false);
    general.set("Scale for Epsilon Active Sets",            1.0  );
    
    /* ----- INEXACT OBJECTIVE AND DERIVATIVES ----------------- */
    general.set("Inexact Objective Function",               false);
    general.set("Inexact Gradient",                         false);
    general.set("Inexact Hessian-Times-A-Vector",           false);
    general.set("Scale for Epsilon Active Sets",            1.0  );
    general.set("Print Verbosity",                          0    );
    
    /* ----- BOUND CONSTRAINED CRITICALITY MEASURE --------------*/
    general.set("Projected Gradient Criticality Measure",   false);
 
    /* ===== SECANT INPUTS ============================================================== */  
    PL &secant = general.sublist("Secant");
      secant.set("Type",                          "Limited-Memory BFGS"                   );
      secant.set("Use as Preconditioner",         false                                   );
      secant.set("Use as Hessian",                false                                   );
      secant.set("Maximum Storage",               10                                      );
      secant.set("Barzilai-Borwein",              1                                       );
      secant.set("User Defined Secant Name",      "Unspecified User Defined Secant Method");
    
    /* ===== KRYLOV INPUTS ========================================================== */
    PL &krylov = general.sublist("Krylov");
      krylov.set("Type",                      "Conjugate Gradients"                   );
      krylov.set("Absolute Tolerance",        1.e-4                                   );
      krylov.set("Relative Tolerance",        1.e-2                                   );
      krylov.set("Iteration Limit",           100                                     );
      krylov.set("User Defined Krylov Name",  "Unspecified User Defined Krylov Method");

  /* ===== STEP SUBLIST ============================================== */    
  PL &step = rol->sublist("Step");

    // Used by OptimizationSolver to select a step 
    step.set("Type","Trust Region");

    PL &linesearch = step.sublist("Line Search");
      linesearch.set("Function Evaluation Limit",                 20   );
      linesearch.set("Sufficient Decrease Tolerance",             1e-4 );
      linesearch.set("Use Previous Step Length as Initial Guess", false);
      linesearch.set("Initial Step Size",                         1.0  );
      linesearch.set("User Defined Initial Step Size",            false);
      linesearch.set("Accept Linesearch Minimizer",               false);
      linesearch.set("Accept Last Alpha",                         false);

      /* ===== DESCENT ALGORITHM SPECIFICATION =============== */
      PL &descent = linesearch.sublist("Descent Method");
        descent.set("Type",               "Quasi-Newton Method");
        descent.set("Nonlinear CG Type",  "Oren-Luenberger    ");

      /* ===== CURVATURE CONDITION SPECIFICATION ============================= */
      PL &curvature = linesearch.sublist("Curvature Condition");
        curvature.set("Type",                         "Strong Wolfe Conditions");
        curvature.set("General Parameter",            0.9                      );
        curvature.set("Generalized Wolfe Parameter",  0.6                      );
   
      /* ===== LINE-SEARCH ALGORITHM SPECIFICATION ======================================== */
      PL &lsmethod = linesearch.sublist("Line-Search Method");
        lsmethod.set("Type",                 "Cubic Interpolation");
        lsmethod.set("Backtracking Rate",    0.5                  );    
        lsmethod.set("Bracketing Tolerance", 1.e-8                );
        lsmethod.set("User Defined Line-Search Name", "Unspecified User Defined Line-Search");

        /* ===== BISECTION METHOD =============== */
        PL &bisection = lsmethod.sublist("Bisection");
        bisection.set("Tolerance",        1.e-10);
        bisection.set("Iteration Limit",  1000  );

        /* ===== BRENT'S METHOD ============ */
        PL &brents = lsmethod.sublist("Brent's");
        brents.set("Tolerance",                    1.e-10);
        brents.set("Iteration Limit",              1000  );
        brents.set("Run Test Upon Initialization", true  );

        /* ===== GOLDEN SECTION =================== */
        PL &golden = lsmethod.sublist("Golden Section");
        golden.set("Tolerance",        1.e-10);
        golden.set("Iteration Limit",  1000  );

        /* ===== PATH-BASED TARGET LEVEL ======================= */
        PL &pathtarg = lsmethod.sublist("Path-Based Target Level");
          pathtarg.set("Target Relaxation Parameter", 1.0);
          pathtarg.set("Upper Bound on Path Length",  1.0);

    /* ===== TRUST REGION ================================================== */
    PL &trustregion = step.sublist("Trust Region");
      trustregion.set("Subproblem Solver",                    "Truncated CG" );
      trustregion.set("Subproblem Model",                     "Kelley-Sachs" );
      trustregion.set("Initial Radius",                       10.0           );
      trustregion.set("Maximum Radius",                       5.e3           );
      trustregion.set("Step Acceptance Threshold",            0.05           );
      trustregion.set("Radius Shrinking Threshold",           0.05           );
      trustregion.set("Radius Growing Threshold",             0.9            );
      trustregion.set("Radius Shrinking Rate (Negative rho)", 0.0625         ); 
      trustregion.set("Radius Shrinking Rate (Positive rho)", 0.25           );
      trustregion.set("Radius Growing Rate",                  2.5            );
      trustregion.set("Sufficient Decrease Parameter",        1.e-2          );
      trustregion.set("Safeguard Size",                       1.e8           );

      /* ===== POST-SMOOTHING SPECIFICATION ============= */ 
      PL &smoothing = trustregion.sublist("Post-Smoothing");
        smoothing.set("Function Evaluation Limit", 20    );
        smoothing.set("Initial Step Size",         1.0   );
        smoothing.set("Tolerance",                 0.9999);
        smoothing.set("Rate",                      0.01  );     
 
      /* ===== COLEMAN-LI MODEL INPUTS ============ */
      PL &coleman = trustregion.sublist("Coleman-Li");
        coleman.set("Maximum Step Back",  0.9999);
        coleman.set("Maximum Step Scale", 1.0   );
        coleman.set("Single Reflection",  true  );
 
      /* ===== CONTROLS FOR INEXACTNESS ======== */
      PL &inexact = trustregion.sublist("Inexact");
   
        /* ===== INEXACT OBJECTIVE VALUE UPDATE ============= */
        PL &value = inexact.sublist("Value");
          value.set("Tolerance Scaling",                 1.e-1);
          value.set("Exponent",                          0.9  );
          value.set("Forcing Sequence Initial Value",    1.0  );
          value.set("Forcing Sequence Update Frequency", 10   );
          value.set("Forcing Sequence Reduction Factor", 0.1  );

        /* ===== INEXACT GRADIENT UPDATE ============ */
        PL &gradient = inexact.sublist("Gradient");      
         gradient.set("Tolerance Scaling",       1.e-1);
         gradient.set("Relative Tolerance",      2.0  );

    /* ===== PRIMAL DUAL ACTIVE SET ==================== */
    PL &activeset = step.sublist("Primal Dual Active Set");
      activeset.set("Dual Scaling",                1.0  );
      activeset.set("Iteration Limit",             10   );
      activeset.set("Relative Step Tolerance",     1.e-8); 
      activeset.set("Relative Gradient Tolerance", 1.e-6);

    /* ===== COMPOSITE STEP ==================== */
    PL &composite = step.sublist("Composite Step");
      composite.set("Output Level",  0);

      /* ===== OPTIMALITY SYSTEM SOLVER ======================== */
      PL &ossolver = composite.sublist("Optimality System Solver");
        ossolver.set("Nominal Relative Tolerance",  1e-8);
        ossolver.set("Fix Tolerance",               true);

      /* ===== TANGENTIAL SUBPROBLEM SOLVER ========================= */
      PL &tansolver = composite.sublist("Tangential Subproblem Solver");
        tansolver.set("Iteration Limit",    20  );
        tansolver.set("Relative Tolerance", 1e-2);

    /* ===== AUGMENTED LAGRANGIAN ======================= */
    PL &auglag = step.sublist("Augmented Lagrangian");

      auglag.set("Use Scaled Augmented Lagrangian",  false);
      auglag.set("Level of Hessian Approximation",       0);

      /* ----- PENALTY PARAMETER UPDATE ----------------------------------- */
      auglag.set("Initial Penalty Parameter",                 1.e1          );
      auglag.set("Penalty Parameter Growth Factor",           100.0         );
      auglag.set("Penalty Parameter Reciprocal Lower Bound",  0.1           );
      auglag.set("Maximum Penalty Parameter",                 1.e8          );

      /* ----- OPTIMALITY TOLERANCE UPDATE ------------------------------ */
      auglag.set("Initial Optimality Tolerance",            1.0           );
      auglag.set("Optimality Tolerance Update Exponent",    0.1           );
      auglag.set("Optimality Tolerance Decrease Exponent",  0.9           );

      /* ----- FEASIBILITY TOLERANCE UPDATE ----------------------------- */
      auglag.set("Initial Feasibility Tolerance",           1.0           );
      auglag.set("Feasibility Tolerance Update Exponent",   0.1           );
      auglag.set("Feasibility Tolerance Decrease Exponent", 0.9           );

      /* ===== SUBPROBLEM SOLVER ======================================== */
      auglag.set("Print Intermediate Optimization History", false         );
      auglag.set("Subproblem Step Type",                    "Trust Region");
      auglag.set("Subproblem Iteration Limit",              1000          );

    /* ===== MOREAU-YOSIDA PENALTY =================== */
    PL &moreau = step.sublist("Moreau-Yosida Penalty");
      moreau.set("Initial Penalty Parameter",       1e2);
      moreau.set("Penalty Parameter Growth Factor", 1.0);

      /* ===== SUBPROBLEM SOLVER =============== */
      PL &mysub = moreau.sublist("Subproblem");
        mysub.set("Optimality Tolerance",  1.e-12);
        mysub.set("Feasibility Tolerance", 1.e-12);
        mysub.set("Print History",         false );
        mysub.set("Iteration Limit",       200   );        

    /* ===== BUNDLE METHOD =================================== */
    PL &bundle = step.sublist("Bundle");

      /* ----- TRUST-REGION RADIUS UPDATE -------------------- */
      bundle.set("Initial Trust-Region Parameter",       1.e1  );
      bundle.set("Maximum Trust-Region Parameter",       1.e8  );
      bundle.set("Tolerance for Trust-Region Parameter", 1.e-4 );

      /* ----- EPSILON SOLUTION STOPPING CONDITION ----------- */      
      bundle.set("Epsilon Solution Tolerance",           1.e-12);

      /* ----- SERIOUS STEP PARAMETERS ----------------------- */
      bundle.set("Upper Threshold for Serious Step",     1.e-1 );
      bundle.set("Lower Threshold for Serious Step",     2.e-1 );
      bundle.set("Upper Threshold for Null Step",        9.e-1 );

      /* ----- BUNDLE INFORMATION ---------------------------- */
      bundle.set("Distance Measure Coefficient",         1.e-6 );
      bundle.set("Maximum Bundle Size",                  50    );
      bundle.set("Removal Size for Bundle Update",       2     );

      /* ----- CUTTING PLANE SUBPROBLEM SOLVER --------------- */
      bundle.set("Cutting Plane Tolerance",              1.e-8 );
      bundle.set("Cutting Plane Iteration Limit",        1000  );

  /* ===== STATUS TEST PARAMETERS ============ */
  PL &status = rol->sublist("Status Test");     
    status.set("Gradient Tolerance",     1.e-10);
    status.set("Constraint Tolerance",   1.e-10);
    status.set("Step Tolerance",         1.e-14);
    status.set("Iteration Limit",        1000  );
 

  return rol;
}


/* SOL Parameters */

inline ROL::Ptr<const ROL::ParameterList> getValidSOLParameters() {
  
  typedef ROL::ParameterList PL;

  ROL::Ptr<PL> sol = ROL::makePtr<PL>("SOL");

  sol->set("Type",  "Risk Neutral");
  sol->set("Store Sampled Value and Gradient", true);
    
  /* ===== RISK MEASURE ============== */
  PL &risk = sol->sublist("Risk Measure");
    risk.set("Name","CVaR");

    /* ===== BPOE ================= */
    PL &bpoe = risk.sublist("bPOE");
      bpoe.set("Moment Order",   2.0);
      bpoe.set("Threshold",      1.0);

    /* ===== EXPONENTIAL UTILITY =============== */
    PL &expo = risk.sublist("Exponential Utility");
      expo.set("Rate", 2.0);

    /* ===== KL DIVERGENCE ================ */
    PL &kldiv = risk.sublist("KL Divergence");
      kldiv.set("Threshold",1.e-2);
  
    /* ===== CHI-SQUARED DIVERGENCE ===== */
    PL &fdiv = risk.sublist("F-Divergence");
      fdiv.set("Threshold",1.e-2);

    /* ===== CVAR SUBLIST =============== */
    PL &cvar = risk.sublist("CVaR");
      cvar.set("Confidence Level",               0.8);
      cvar.set("Convex Combination Parameter",   0.8);
      cvar.set("Smoothing Parameter",          1.e-2);

      PL &cvar_dist = cvar.sublist("Distribution");
  
        cvar_dist.set("Name",  "Parabolic");

        PL &cvar_para = cvar_dist.sublist("Parbolic");
          cvar_para.set("Lower Bound",-0.5);
          cvar_para.set("Upper Bound", 0.5);
    
    /* ===== HMCR SUBLIST =============== */
    PL &hmcr = risk.sublist("HMCR");
      hmcr.set("Confidence Level",              0.8  );
      hmcr.set("Convex Combination Parameter",  0.8  );
      hmcr.set("Order",                         2    );
      hmcr.set("Smoothing Parameter",           1.e-2);

      PL &hmcr_dist = hmcr.sublist("Distribution");

        hmcr_dist.set("Name", "Dirac");

        PL &hmcr_dirac = hmcr_dist.sublist("Dirac");
          hmcr_dirac.set("Location",0.0);



  return sol;
}


} // namespace ROL

#endif // ROL_VALID_PARAMETERS_H



