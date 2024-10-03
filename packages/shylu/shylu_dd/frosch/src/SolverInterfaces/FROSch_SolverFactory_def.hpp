// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SOLVER_FACTORY_DEF_HPP
#define _FROSCH_SOLVER_FACTORY_DEF_HPP

#include "FROSch_SolverFactory_decl.hpp"

#include "Stratimikos_FROSch_def.hpp"


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC, class LO, class GO, class NO>
    typename SolverFactory<SC,LO,GO,NO>::SolverPtr SolverFactory<SC,LO,GO,NO>::Build(ConstXMatrixPtr k,
                                                                                     ParameterListPtr parameterList,
                                                                                     string description)
    {
        const string solverType = parameterList->get("SolverType","Amesos2");
        if (!solverType.compare("Amesos")) {
#ifdef HAVE_SHYLU_DDFROSCH_AMESOS
            FROSCH_ASSERT(k->getRowMap()->lib()==UseEpetra,"FROSch::SolverFactory: Amesos is not compatible with Tpetra.");
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
            return AmesosSolverEpetraPtr(new AmesosSolverEpetra<SC,LO,GO,NO>(k,parameterList,description));
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Epetra");
#endif
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Amesos");
#endif
        } else if (!solverType.compare("Amesos2")) {
            if (k->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                return Amesos2SolverEpetraPtr(new Amesos2SolverEpetra<SC,LO,GO,NO>(k,parameterList,description));
#else
                ThrowErrorMissingPackage("FROSch::SolverFactory","Epetra");
#endif
            } else if (k->getRowMap()->lib()==UseTpetra) {
                return Amesos2SolverTpetraPtr(new Amesos2SolverTpetra<SC,LO,GO,NO>(k,parameterList,description));
            } else {
                FROSCH_ASSERT(false, "FROSch::SolverFactory: This can't happen. Either use Epetra or Tetra linear algebra stack.");
            }
        } else if (!solverType.compare("Belos")) {
#ifdef HAVE_SHYLU_DDFROSCH_BELOS
            if (k->getRowMap()->lib()==UseEpetra) {
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
                return BelosSolverEpetraPtr(new BelosSolverEpetra<SC,LO,GO,NO>(k,parameterList,description));
#else
                ThrowErrorMissingPackage("FROSch::SolverFactory","Epetra");
#endif
            } else if (k->getRowMap()->lib()==UseTpetra) {
                return BelosSolverTpetraPtr(new BelosSolverTpetra<SC,LO,GO,NO>(k,parameterList,description));
            } else {
                FROSCH_ASSERT(false, "FROSch::SolverFactory: This can't happen. Either use Epetra or Tetra linear algebra stack.");
            }
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Belos");
#endif
        } else if (!solverType.compare("FROSchPreconditioner")) {
            return FROSchPreconditionerPtr(new FROSchPreconditioner<SC,LO,GO,NO>(k,parameterList,description));
        } else if (!solverType.compare("Ifpack2")) {
#ifdef HAVE_SHYLU_DDFROSCH_IFPACK2
            FROSCH_ASSERT(k->getRowMap()->lib()==UseTpetra,"FROSch::SolverFactory: Ifpack2 is not compatible with Epetra.");
            return Ifpack2PreconditionerTpetraPtr(new Ifpack2PreconditionerTpetra<SC,LO,GO,NO>(k,parameterList,description));
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Ifpack2");
#endif
        } else if (!solverType.compare("MueLu")) {
#ifdef HAVE_SHYLU_DDFROSCH_MUELU
            return MueLuPreconditionerPtr(new MueLuPreconditioner<SC,LO,GO,NO>(k,parameterList,description));
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","MueLu");
#endif
        } else if (!solverType.compare("ThyraPreconditioner")) {
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
#ifdef HAVE_SHYLU_DDFROSCH_STRATIMIKOS
            return ThyraPreconditionerPtr(new ThyraPreconditioner<SC,LO,GO,NO>(k,parameterList,description));
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Stratimikos");
#endif
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Thyra");
#endif
        } else if (!solverType.compare("ThyraSolver")) {
#ifdef HAVE_SHYLU_DDFROSCH_THYRA
#ifdef HAVE_SHYLU_DDFROSCH_STRATIMIKOS
            return ThyraSolverPtr(new ThyraSolver<SC,LO,GO,NO>(k,parameterList,description));
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Stratimikos");
#endif
#else
            ThrowErrorMissingPackage("FROSch::SolverFactory","Thyra");
#endif
        } else {
            FROSCH_ASSERT(false,"SolverType unknown...");
        }
        return null; // Can we get rid of this? Without there used to be a warning...
    }
}

#endif
