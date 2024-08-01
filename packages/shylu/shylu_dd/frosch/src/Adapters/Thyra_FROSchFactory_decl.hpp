// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _THYRA_FROSCH_FACTORY_DECL_HPP
#define _THYRA_FROSCH_FACTORY_DECL_HPP

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

//Tpetra
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
#include <Xpetra_TpetraHalfPrecisionOperator.hpp>
#endif

//Epetra
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_config.h>
#endif

//FROSch
#include <Thyra_FROSchLinearOp_def.hpp>
#include <FROSch_Tools_def.hpp>

#include <FROSch_SchwarzPreconditioners_fwd.hpp>


namespace Thyra {

    using namespace FROSch;
    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,
              class LO,
              class GO,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSchFactory : public Thyra::PreconditionerFactoryBase<SC> {

    protected:

        using CommPtr                       = RCP<const Comm<int> >;

        using LinearOpBasePtr               = RCP<LinearOpBase<SC> >;
        using ConstLinearOpBasePtr          = RCP<const LinearOpBase<SC> >;
        using ConstLinearOpSourceBasePtr    = RCP<const LinearOpSourceBase<SC> >;

        using ConstVectorSpaceBasePtr       = RCP<const VectorSpaceBase<SC> >;

        using PreconditionerBasePtr         = RCP<PreconditionerBase<SC> >;

        using XMap                          = Map<LO,GO,NO>;
        using ConstXMap                     = const XMap;
        using XMapPtr                       = RCP<XMap>;
        using ConstXMapPtr                  = RCP<ConstXMap>;
        using XMapPtrVecPtr                 = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr            = ArrayRCP<ConstXMapPtr>;

        using XMatrix                       = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                    = RCP<XMatrix>;
        using ConstXMatrixPtr               = RCP<const XMatrix>;

        using XCrsMatrix                    = CrsMatrix<SC,LO,GO,NO>;
        using XCrsMatrixPtr                 = RCP<XCrsMatrix>;
        using ConstXCrsMatrixPtr            = RCP<const XCrsMatrix>;

        using XMultiVector                  = MultiVector<SC,LO,GO,NO>;
        using ConstXMultiVector             = const XMultiVector;
        using XMultiVectorPtr               = RCP<XMultiVector>;
        using ConstXMultiVectorPtr          = RCP<ConstXMultiVector>;
        using XMultiVectorPtrVecPtr         = ArrayRCP<XMultiVectorPtr>;
        using ConstXMultiVectorPtrVecPtr    = ArrayRCP<ConstXMultiVectorPtr>;

#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
        using HalfPrecOp                    = Xpetra::TpetraHalfPrecisionOperator<SC,LO,GO,NO>;
        using HalfSC                        = typename HalfPrecOp::HalfScalar;
        using HalfPrecMatrix                = Matrix<HalfSC,LO,GO,NO>;
        using HalfPrecMultiVector           = MultiVector<HalfSC,LO,GO,NO>;
#endif

        using ParameterListPtr              = RCP<ParameterList>;
        using ConstParameterListPtr         = RCP<const ParameterList>;

        using DofOrderingVecPtr             = ArrayRCP<DofOrdering>;

        using UN                            = unsigned;

        using GOVecPtr                      = ArrayRCP<GO> ;

        using SCVecPtr                      = ArrayRCP<SC>;

        using UNVecPtr                      = ArrayRCP<UN>;

        using LOVecPtr                      = ArrayRCP<LO>;

    public:

        //Constructor
        FROSchFactory();

        //Overridden from PreconditionerFactory Base
        bool isCompatible(const LinearOpSourceBase<SC>& fwdOp) const;

        PreconditionerBasePtr createPrec() const;

        void initializePrec(const ConstLinearOpSourceBasePtr& fwdOpSrc,
                            PreconditionerBase<SC>* prec,
                            const ESupportSolveUse supportSolveUse) const;

        void uninitializePrec(PreconditionerBase<SC>* prec,
                              ConstLinearOpSourceBasePtr* fwdOp,
                              ESupportSolveUse* supportSolveUse) const;

        void setParameterList(const ParameterListPtr& paramList);

        ParameterListPtr unsetParameterList();

        ParameterListPtr getNonconstParameterList();

        ConstParameterListPtr getParameterList() const;

        ConstParameterListPtr getValidParameters() const;

        string description() const;

    private:

        ConstXMapPtr extractRepeatedMap(CommPtr comm,
                                        UnderlyingLib lib) const;

        ConstXMultiVectorPtr extractCoordinatesList(CommPtr comm,
                                                    UnderlyingLib lib) const;

        ConstXMultiVectorPtr extractNullSpace(CommPtr comm,
                                              UnderlyingLib lib) const;


        ParameterListPtr paramList_ = rcp(new ParameterList());
    };

}

#endif
