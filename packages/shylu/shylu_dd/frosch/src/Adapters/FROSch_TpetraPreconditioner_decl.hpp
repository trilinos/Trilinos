// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TPETRAPRECONDITIONER_DECL_HPP
#define _FROSCH_TPETRAPRECONDITIONER_DECL_HPP

#if defined(HAVE_XPETRA_TPETRA)
#include <Tpetra_Operator.hpp>

#include <Xpetra_Operator.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_Types.h>
#include <FROSch_SchwarzPreconditioner_decl.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class TpetraPreconditioner : public Tpetra::Operator<SC,LO,GO,NO> {

    protected:
        using TMap = Tpetra::Map<LO,GO,NO>;
        using TMultiVector = Tpetra::MultiVector<SC,LO,GO,NO>;
    public:

        TpetraPreconditioner();
        TpetraPreconditioner(Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> > preconditioner);

        ~TpetraPreconditioner();

        int initialize(bool useDefaultParameters = true);

        int compute();
        
        // Y = alpha * A^mode * X + beta * Y
        void apply(const TMultiVector &X,
                         TMultiVector &Y,
                   ETransp mode=NO_TRANS,
                   SC alpha=ScalarTraits<SC>::one(),
                   SC beta=ScalarTraits<SC>::zero()) const;

        Teuchos::RCP<const TMap> getDomainMap() const;

        Teuchos::RCP<const TMap> getRangeMap() const;

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const;

        bool isInitialized() const;

        bool isComputed() const;

        void residual(const TMultiVector & X,
                      const TMultiVector & B,
                            TMultiVector& R) const;

        Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> >
        getSchwarzPreconditioner();

    protected:
        Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> > preconditioner_;
    };

}
#endif

#endif
