//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef SHYLU_CONFIG_H
#define SHYLU_CONFIG_H

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
    string libName;             // library for the outer solver
    string schurSolver;         // Solver for the Schur complement
    string schurAmesosSolver;   // Amesos solver for the Schur complement
    string diagonalBlockSolver; // Solver to use to factorize the diagonal blocks
    string schurPreconditioner; // Preconditioner for the inner iterations on Sbar (AztecOO-Exact)
    bool silent_subiter;
    int sep_type;
    int debug_level;
    //DebugManager dm;
    int reset_iter;             // When should we reset the guided_probing
    int overlap;
    bool amesosForDiagonal;
} shylu_config;

#endif // SHYLU_CONFIG_H
