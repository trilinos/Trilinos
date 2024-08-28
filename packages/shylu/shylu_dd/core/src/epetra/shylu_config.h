// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLU_CONFIG_H
#define SHYLU_CONFIG_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

typedef struct
{
    int sym;                    // flag for symmetry
    double Sdiagfactor;         // % of diagonals added to Schur complement
    int schurApproxMethod;      // ==1 implies blockdiagonal + A22
                                // ==2 implies dropping based
    							// ==4 implies IQR
    double relative_threshold;  // Relative threshold for dropping
                                // only used if schurApproxMethod == 2
    int inner_maxiters;         // maximum iterations for inner solver
    double inner_tolerance;     // relative residual tolerance for inner solver
    std::string libName;             // library for the outer solver
    std::string schurSolver;         // Solver for the Schur complement
    std::string schurAmesosSolver;   // Amesos solver for the Schur complement
    std::string diagonalBlockSolver; // Solver to use to factorize the diagonal blocks
    std::string schurPreconditioner; // Preconditioner for the inner iterations on Sbar (AztecOO-Exact)
    bool silent_subiter;
    int sep_type;
    int debug_level;
    //DebugManager dm;
    int reset_iter;             // When should we reset the guided_probing
    int overlap;
    bool amesosForDiagonal;
} shylu_config;

#endif // SHYLU_CONFIG_H
