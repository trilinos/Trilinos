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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_SOLVER_FACTORY_DECL_HPP
#define _FROSCH_SOLVER_FACTORY_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

// FROSch
#include <FROSch_Solver_def.hpp>
#if defined(HAVE_SHYLU_DDFROSCH_AMESOS) && defined(HAVE_SHYLU_DDFROSCH_EPETRA)
#include <FROSch_AmesosSolverEpetra_def.hpp>
#endif
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <FROSch_Amesos2SolverEpetra_def.hpp>
#endif
#include <FROSch_Amesos2SolverTpetra_def.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <FROSch_BelosSolverEpetra_def.hpp>
#endif
#include <FROSch_BelosSolverTpetra_def.hpp>
#endif
#include <FROSch_FROSchPreconditioner_def.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
#include <FROSch_Ifpack2PreconditionerTpetra_def.hpp>
#endif
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
#include <FROSch_MueLuPreconditioner_def.hpp>
#endif
#if defined(HAVE_SHYLU_DDFROSCH_THYRA) && defined(HAVE_SHYLU_DDFROSCH_STRATIMIKOS)
#include <FROSch_ThyraPreconditioner_def.hpp>
#include <FROSch_ThyraSolver_def.hpp>
#endif


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class SolverFactory {

    protected:

        using XMatrixPtr                        = typename Solver<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename Solver<SC,LO,GO,NO>::ConstXMatrixPtr;

        using ParameterListPtr                  = typename Solver<SC,LO,GO,NO>::ParameterListPtr;

        using SolverPtr                         = RCP<Solver<SC,LO,GO,NO> >;
#if defined(HAVE_SHYLU_DDFROSCH_AMESOS) && defined(HAVE_SHYLU_DDFROSCH_EPETRA)
        using AmesosSolverEpetraPtr             = RCP<AmesosSolverEpetra<SC,LO,GO,NO> >;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
        using Amesos2SolverEpetraPtr            = RCP<Amesos2SolverEpetra<SC,LO,GO,NO> >;
#endif
        using Amesos2SolverTpetraPtr            = RCP<Amesos2SolverTpetra<SC,LO,GO,NO> >;
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
        using BelosSolverEpetraPtr              = RCP<BelosSolverEpetra<SC,LO,GO,NO> >;
#endif
        using BelosSolverTpetraPtr              = RCP<BelosSolverTpetra<SC,LO,GO,NO> >;
#endif
        using FROSchPreconditionerPtr           = RCP<FROSchPreconditioner<SC,LO,GO,NO> >;
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
        using Ifpack2PreconditionerTpetraPtr    = RCP<Ifpack2PreconditionerTpetra<SC,LO,GO,NO> >;
#endif
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
        using MueLuPreconditionerPtr            = RCP<MueLuPreconditioner<SC,LO,GO,NO> >;
#endif
#if defined(HAVE_SHYLU_DDFROSCH_THYRA) && defined(HAVE_SHYLU_DDFROSCH_STRATIMIKOS)
        using ThyraPreconditionerPtr            = RCP<ThyraPreconditioner<SC,LO,GO,NO> >;
        using ThyraSolverPtr                    = RCP<ThyraSolver<SC,LO,GO,NO> >;
#endif

    public:

        static SolverPtr Build(ConstXMatrixPtr k,
                               ParameterListPtr parameterList,
                               string description);
    };

}

#endif
