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

#ifndef THYRA_FROSCH_XPETRA_FACTORY_DECL_HPP
#define THYRA_FROSCH_XPETRA_FACTORY_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

//Thyra
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include "Thyra_EpetraLinearOp.hpp"
#endif
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

//Teuchos
#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Array.hpp"
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

//Xpetra
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_ThyraUtils.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>
#endif

//Epetra
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_config.h>
#endif

//FROSch
#include <FROSch_AlgebraicOverlappingPreconditioner_def.hpp>
#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_RGDSWPreconditioner_def.hpp>
#include <FROSch_OneLevelPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
#include <FROSch_TwoLevelBlockPreconditioner_def.hpp>
#include <Thyra_FROSchLinearOp_def.hpp>
#include <FROSch_Tools_def.hpp>

#include "Kokkos_DefaultNode.hpp"

namespace Thyra {
    
    using namespace FROSch;
    using namespace Teuchos;
    
    template <class SC,class LO,class GO,class NO=KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSchFactory : public Thyra::PreconditionerFactoryBase<SC> {
        
    public:

        typedef Teuchos::ArrayRCP<DofOrdering> DofOrderingVecPtr;
        
        typedef unsigned UN;
        
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;
        
        typedef Teuchos::ArrayRCP<UN> UNVecPtr;
        
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        
        // More typedefs!!!
        
        //Constructor
        FROSchFactory();
        
        //Overridden from PreconditionerFactory Base
        bool isCompatible(const LinearOpSourceBase<SC>& fwdOp) const;
        
        Teuchos::RCP<PreconditionerBase<SC> > createPrec() const;
        
        void initializePrec(const RCP<const LinearOpSourceBase<SC> >& fwdOpSrc,
                            PreconditionerBase<SC>* prec,
                            const ESupportSolveUse supportSolveUse) const;
        
        void uninitializePrec(PreconditionerBase<SC>* prec,
                              RCP<const LinearOpSourceBase<SC> >* fwdOp,
                              ESupportSolveUse* supportSolveUse) const;
        
        void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
        
        RCP<ParameterList> unsetParameterList();
        
        RCP<ParameterList>getNonconstParameterList();
        
        RCP<const ParameterList> getParameterList() const;
        
        RCP<const ParameterList> getValidParameters() const;
        
        std::string description() const;
        private:
        Teuchos::RCP<ParameterList> paramList_;
                
    };
}

    

#endif


