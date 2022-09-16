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
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class TpetraPreconditioner : public Tpetra::Operator<SC,LO,GO,NO> {

    protected:
        using TMap = Tpetra::Map<LO,GO,NO>;
        using TMultiVector = Tpetra::MultiVector<SC,LO,GO,NO>;
    public:

        TpetraPreconditioner();
        TpetraPreconditioner(Teuchos::RCP<SchwarzPreconditioner<SC,LO,GO,NO> > preconditioner);

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

        Teuchos::RCP<SchwarzPreconditioner<SC,LO,GO,NO> >
        getSchwarzPreconditioner();

    protected:
        Teuchos::RCP<SchwarzPreconditioner<SC,LO,GO,NO> > preconditioner_;
    };

}
#endif

#endif
