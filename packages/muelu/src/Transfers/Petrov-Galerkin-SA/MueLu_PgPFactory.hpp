/*
 * MueLu_PgPFactory.hpp
 *
 *  Created on: 23.09.2011
 *      Author: tobias
 */

#ifndef MUELU_PGPFACTORY_HPP_
#define MUELU_PGPFACTORY_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_TpetraVector.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
    @class PgPFactory class.
    @brief Factory for building Petrov-Galerkin Smoothed Aggregation prolongators.
 */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class PgPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, PgPFactory<AA,BB,CC,DD, EE> &factory);

public:

  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
   */
  PgPFactory(RCP<PFactory> InitialPFact = Teuchos::null, RCP<SingleLevelFactoryBase> AFact = Teuchos::null)
  : initialPFact_(InitialPFact), AFact_(AFact),
    diagonalView_("current") {
    if(initialPFact_ == Teuchos::null)
    {
      // use tentative P factory as default
      initialPFact_ = rcp(new TentativePFactory());
    }
  }

  //! Destructor.
  virtual ~PgPFactory() {}

  //@}

  //! @name Set methods.
  //@{

  //! Change view of diagonal.
  void SetDiagonalView(std::string const& diagView) {
    diagonalView_ = diagView;
  }
  //@}

  //! @name Get methods.
  //@{

  //! Returns current view of diagonal.
  std::string GetDiagonalView() {
    return diagonalView_;
  }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    //coarseLevel.Request("P",initialPFact_.get());
    //fineLevel.Request("A",AFact_.get());
    fineLevel.DeclareInput("A",AFact_.get());
    coarseLevel.DeclareInput("P",initialPFact_.get());
  };
  //@}

  //! @name Build methods.
  //@{

  /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
   */
  void Build(Level& fineLevel, Level &coarseLevel) const {

    std::ostringstream buf; buf << coarseLevel.GetLevelID();
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP_"+buf.str()));
    timer->start(true);

    // Level Get
    RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_.get());
    RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_.get());

    if(restrictionMode_)
      A = Utils2<Scalar,LocalOrdinal,GlobalOrdinal>::Transpose(A,true); // build transpose of A explicitely

    //fineLevel.Release("A", AFact_.get());
    //coarseLevel.Release("P", initialPFact_.get());

    std::cout << "after Ptent has been calculated" << std::endl;

    // calculate D^{-1} * A * Ptent

    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> DinvAP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

    doFillComplete=true;
    optimizeStorage=false;
    Teuchos::ArrayRCP<Scalar> diag = Utils::GetMatrixDiagonal(A);
    Utils::MyOldScaleMatrix(DinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    // calculate local damping factors

    // compute D^{-1} * A * D^{-1} * A * P0
    doFillComplete=true;
    optimizeStorage=false;
    RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAP0,false,doFillComplete,optimizeStorage);

    doFillComplete=true;
    optimizeStorage=false;
    Utils::MyOldScaleMatrix(DinvADinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > colbasedomegamap = BuildLocalReplicatedColMap(DinvAP0);

    std::cout<< "after BuildLocallyReplicatedColMap" << std::endl;

    //if(colbasedomegamap->isDistributed() == true) std::cout << "colbasedomegamap is distributed" << std::endl; //throw("colbasedomegamap is distributed globally?");
    if(colbasedomegamap->getMinAllGlobalIndex() != DinvAP0->getDomainMap()->getMinAllGlobalIndex()) std::cout << "MinAllGID does not match" << std::endl; //throw("MinAllGID does not match");
    if(colbasedomegamap->getMaxAllGlobalIndex() != DinvAP0->getDomainMap()->getMaxAllGlobalIndex()) std::cout << "MaxAllGID does not match" << std::endl; //throw("MaxAllGID does not match");
    //if(colbasedomegamap->getGlobalNumElements() != DinvAP0->getDomainMap()->getGlobalNumElements()) std::cout << "NumGlobalElements do not match" << std::endl; //throw("NumGlobalElements do not match");

    RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Numerator = Teuchos::null;
        //Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colbasedomegamap);
    RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Denominator = Teuchos::null;
        //Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colbasedomegamap);
    RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ColBasedOmegas = Teuchos::null;
       // Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colbasedomegamap);

    std::cout<< "do MultiplyAll" << std::endl;
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

    Numerator = MultiplyAll(DinvAP0, DinvADinvAP0, colbasedomegamap);
    Denominator = MultiplySelfAll(DinvADinvAP0, colbasedomegamap);

    // check for zeros in denominator -> error
    size_t zeros_in_denominator = 0;
    Teuchos::ArrayRCP< const Scalar > Numerator_data   = Numerator->getData(0);
    Teuchos::ArrayRCP< const Scalar > Denominator_data = Denominator->getData(0);
    for(size_t i=0; i<Denominator_data.size(); i++)
    {
      if(Denominator_data[i] == Teuchos::ScalarTraits<Scalar>::zero()) zeros_in_denominator ++;
    }
    if(zeros_in_denominator>Teuchos::ScalarTraits<Scalar>::zero()) std::cout << "There are " << zeros_in_denominator<< " zeros in Denominator. very suspicious!" << std::endl;

    // build ColBasedOmegas
    if (DinvAP0->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      RCP<Xpetra::EpetraVector> eColBasedOmegas =
                rcp(new Xpetra::EpetraVector(colbasedomegamap,true));

      for(size_t i=0; i<Denominator_data.size(); i++)
      {
        eColBasedOmegas->getEpetra_Vector()->ReplaceMyValue(i,0,Numerator_data[i]/Denominator_data[i]);
      }
      ColBasedOmegas = eColBasedOmegas;

#else
throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }

    // Tpetra-specific code
    else if(DinvAP0->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tColBasedOmegas =
                rcp(new Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(colbasedomegamap,true));
      for(size_t i=0; i<Denominator_data.size(); i++)
      {
        tColBasedOmegas->getTpetra_Vector()->replaceLocalValue(i, Numerator_data[i]/Denominator_data[i]);
      }
      ColBasedOmegas = tColBasedOmegas;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }

    // check for negative entries in ColBasedOmegas
    Teuchos::ArrayRCP< Scalar > ColBasedOmega_data   = ColBasedOmegas->getDataNonConst(0);
    for(size_t i=0; i<ColBasedOmega_data.size(); i++)
    {
      if(ColBasedOmega_data[i] < Teuchos::ScalarTraits<Scalar>::zero()) ColBasedOmega_data[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }

    // create RowBasedOmegas
    RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RowBasedOmegas = Teuchos::null;

    // build ColBasedOmegas
    if (DinvAP0->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      RCP<Xpetra::EpetraVector> eRowBasedOmegas =
                rcp(new Xpetra::EpetraVector(DinvAP0->getRowMap(),true));
      eRowBasedOmegas->putScalar(-666.0);

      for(size_t row=0; row<DinvAP0->getNodeNumRows(); row++)
      {
        size_t nnz = DinvAP0->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const GlobalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        DinvAP0->getLocalRowView(row, indices, vals);

        if(indices.size() != nnz) std::cout << "number of nonzeros not equal to number of indices?" << std::endl;

        if(nnz == 0)
        {
          eRowBasedOmegas->getEpetra_Vector()->ReplaceMyValue(row,0,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else if(nnz == 1)
        {
          // dirichlet dofs
          eRowBasedOmegas->getEpetra_Vector()->ReplaceMyValue(row,0,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else
        {
          for(size_t j=0; j<nnz; j++)
          {
            GlobalOrdinal col_gid = indices[j]; // local id = global id
            LocalOrdinal col_omega_lid = colbasedomegamap->getLocalElement(col_gid);
            Scalar omega = ColBasedOmega_data[col_omega_lid];
            Epetra_Vector* eeRowBasedOmegas = eRowBasedOmegas->getEpetra_Vector();
            if((*eeRowBasedOmegas)[row] == -666.0)
              eeRowBasedOmegas->ReplaceMyValue(row,0,omega);
            else if (omega < (*eeRowBasedOmegas)[row])
              eeRowBasedOmegas->ReplaceMyValue(row,0,omega);
          }
        }
      }
      RowBasedOmegas = eRowBasedOmegas;

#else
throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }

    // Tpetra-specific code
    else if(DinvAP0->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tRowBasedOmegas =
                rcp(new Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(DinvAP0->getRowMap(),true));
      tRowBasedOmegas->putScalar(-666.0);
      for(size_t row=0; row<DinvAP0->getNodeNumRows(); row++)
      {
        size_t nnz = DinvAP0->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const GlobalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        DinvAP0->getLocalRowView(row, indices, vals);

        if(indices.size() != nnz) std::cout << "number of nonzeros not equal to number of indices?" << std::endl;

        if(nnz == 0)
        {
          tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else if(nnz == 1)
        {
          // dirichlet dofs
          tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else
        {
          for(size_t j=0; j<nnz; j++)
          {
            GlobalOrdinal col_gid = indices[j]; // local id = global id
            LocalOrdinal col_omega_lid = colbasedomegamap->getLocalElement(col_gid);
            Scalar omega = ColBasedOmega_data[col_omega_lid];
            RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ttRowBasedOmegas = tRowBasedOmegas->getTpetra_Vector();
            Teuchos::ArrayRCP<const Scalar> localRowData = ttRowBasedOmegas->getData(0);
            if(localRowData[row] == -666.0)
              tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,omega);
            else if (omega < localRowData[row])
              tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,omega);
          }
        }
      }

      RowBasedOmegas = tRowBasedOmegas;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }

    // check for negative entries in RowBasedOmegas
    Teuchos::ArrayRCP< Scalar > RowBasedOmega_data = RowBasedOmegas->getDataNonConst(Teuchos::as<size_t>(0));
    for(size_t i=0; i<RowBasedOmega_data.size(); i++)
    {
      if(RowBasedOmega_data[i] < Teuchos::ScalarTraits<Scalar>::zero()) RowBasedOmega_data[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }


    ///////////////////



    doFillComplete=true;
    optimizeStorage=false;
    Utils::MyOldScaleMatrix(DinvAP0,RowBasedOmega_data,false,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    RCP<Operator> P_smoothed = Teuchos::null;
    Utils::TwoMatrixAdd(Ptent, false, Teuchos::ScalarTraits<Scalar>::one(),
                                         DinvAP0, false, Teuchos::ScalarTraits<Scalar>::one(),
                                         P_smoothed);
    P_smoothed->fillComplete(Ptent->getDomainMap(), Ptent->getRangeMap());

    P_smoothed->describe(*fos,Teuchos::VERB_EXTREME);

    ////////////////////

    // Level Set
    if(!restrictionMode_)
    {
        // prolongation factory is in prolongation mode
        coarseLevel.Set("P", P_smoothed, this);
    }
    else
    {
        // prolongation factory is in restriction mode
        RCP<Operator> R = Utils2<SC,LO,GO>::Transpose(P_smoothed,true); // use Utils2 -> specialization for double
        coarseLevel.Set("R", R, this);
    }

    //Test(DinvAP0);
    //rcolVec->putScalar(Teuchos::as<const Scalar>(2.0));
    //localomegaVec->doImport(rcolVec,omegaImport);
    //localomegaVec->describe(*fos,Teuchos::VERB_EXTREME);
    /////////////////////

    /*Teuchos::ArrayView<const LocalOrdinal> lindices_right;
        Teuchos::ArrayView<const Scalar> lvals_right;

        DinvADinvAP0->getLocalRowView(49, lindices_right, lvals_right);*/
    //std::cout << "localindices for row l 49 on proc " << comm->getRank() << " " <<  lindices_right << std::endl;
    //for(size_t i = 0; i< lindices_right.size(); i++)
    //    std::cout << "PROC " << comm->getRank() << " LID: " << loclist[lindices_right[i]] << " GID: " << DinvADinvAP0->getColMap()->getGlobalElement(loclist[lindices_right[i]]) << std::endl;
    //std::cout << "localvals for row l 49 on proc " << comm->getRank() << " " <<  lvals_right << std::endl;

    ///////////////// minimize with respect to the (D^{-1} A)' D^{-1} A norm.
    //
    //               diag( R0 (A D^{-1}' D^{-1} A' D^{-1} A' R0' )
    //  omega = ---------------------------------------------------------
    //           diag( R0 A D^{-1} A D^{-1} D^{-1} A' D^{-1} A' R0' )
    //Teuchos::Array<Scalar> Numerator(localreplicatedcolgids.size(),0.0);
    //MultiplyAll(DinvAP0,DinvADinvAP0,Numerator);  // -> DinvAP0_subset

    //std::cout << Numerator << std::endl;
    // smooth Pgi

    // do wiggle-cols handling?

  }

  void BuildP(Level &fineLevel, Level &coarseLevel) const {
    std::cout << "TODO: remove me" << std::endl;
  }

  //@}


  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MultiplySelfAll(const RCP<Operator>& Op, const RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >& InnerProdMap) const
  {
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

    // test
    /*Teuchos::RCP<Xpetra::Operator<double> > OpT = MueLu::Utils2<double>::Transpose(Op,false);
    Teuchos::RCP<Xpetra::Operator<double> > OpTOp = MueLu::Utils<double>::TwoMatrixMultiply(OpT,false,Op,false);
    OpTOp->describe(*fos,Teuchos::VERB_EXTREME);*/


    // Epetra-specific code
    if (Op->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT


      Epetra_Map eInnerProdMap = toEpetra(InnerProdMap);
      RCP<Epetra_Vector> InnerProd_local = rcp(new Epetra_Vector(eInnerProdMap,true));

      for(size_t n=0; n<Op->getNodeNumRows(); n++)
      {
        Teuchos::ArrayView<const LocalOrdinal> lindices;
        Teuchos::ArrayView<const Scalar> lvals;
        Op->getLocalRowView(n, lindices, lvals);

        std::vector<Scalar> retvals;
        std::vector<GlobalOrdinal> retindx;
        retvals.reserve(lindices.size());  // we have not more than nnz_left possible matchings of left and right gids
        retindx.reserve(lindices.size());

        for(size_t i=0; i<Teuchos::as<size_t>(lindices.size()); i++)
        {
          GlobalOrdinal gid = Op->getColMap()->getGlobalElement(lindices[i]);
          retvals.push_back(lvals[i]*lvals[i]);
          retindx.push_back(gid); // note: we save the gids
        }

        InnerProd_local->SumIntoGlobalValues(Teuchos::as<size_t>(retvals.size()),&retvals[0],&retindx[0]); // we use the gids
      }

      // safety check
      if(InnerProd_local->GlobalLength() != InnerProd_local->MyLength()) std::cout << "Global and Local length have to be the same!" << std::endl;

      /////////////// extract all local gids and values
      std::vector<Scalar>           localvalues(InnerProd_local->GlobalLength(),Teuchos::ScalarTraits<Scalar>::zero());
      std::vector<GlobalOrdinal>    globalindices(InnerProd_local->GlobalLength(),-Teuchos::ScalarTraits<Scalar>::one());
      std::vector<Scalar>           globalvalues(InnerProd_local->GlobalLength(),Teuchos::ScalarTraits<Scalar>::zero());
      for(size_t i=0; i<localvalues.size(); i++)
      {
        localvalues[i] = (*InnerProd_local)[i];
        globalindices[i] = InnerProd_local->Map().GID(i); //i; // store indices (attention: local index = global index)
      }

      /////////////// sum up all values to global
      eInnerProdMap.Comm().SumAll(&localvalues[0],&globalvalues[0],localvalues.size()); // sum up all local vectors

      /////////////// save all global information in local replicated vector InnerProd_local
      InnerProd_local->ReplaceGlobalValues(globalvalues.size(),&globalvalues[0],&globalindices[0]);

      /////////////// free some memory
      localvalues.clear();
      globalindices.clear();
      globalvalues.clear();

      RCP<Xpetra::EpetraVector> xInnerProd_local = Teuchos::rcp(new Xpetra::EpetraVector(InnerProd_local));
      return xInnerProd_local;
      #else
      throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }

    // Tpetra-specific code
    else if(Op->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      // Tpetra specific!
      RCP<Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > InnerProd_local =
          //Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(InnerProdMap, true);
          rcp(new Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(InnerProdMap,true));

      for(size_t n=0; n<Op->getNodeNumRows(); n++)
      {
        Teuchos::ArrayView<const LocalOrdinal> lindices;
        Teuchos::ArrayView<const Scalar> lvals;

        Op->getLocalRowView(n, lindices, lvals);

        for(size_t i=0; i<Teuchos::as<size_t>(lindices.size()); i++)
        {
            GlobalOrdinal gid = Op->getColMap()->getGlobalElement(lindices[i]);
            InnerProd_local->sumIntoLocalValue(gid, 0, lvals[i] * lvals[i]);
        }
      }

      InnerProd_local->reduce();

      return InnerProd_local;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }

    return Teuchos::null;
  }

  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MultiplyAll(const RCP<Operator>& left, const RCP<Operator>& right, const RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >& InnerProdMap) const
  {
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

    if(!left->getDomainMap()->isSameAs(*right->getDomainMap()))
      std::cout << "domain maps of left and right do not match" << std::endl;
    if(!left->getRowMap()->isSameAs(*right->getRowMap()))
      std::cout << "row maps of left and right do not match" << std::endl;

    // test
    /*Teuchos::RCP<Xpetra::Operator<double> > leftT = MueLu::Utils2<double>::Transpose(left,false);
    Teuchos::RCP<Xpetra::Operator<double> > leftTright = MueLu::Utils<double>::TwoMatrixMultiply(leftT,false,right,false);
    leftTright->describe(*fos,Teuchos::VERB_EXTREME);*/


    // Epetra-specific code
    if (left->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT


      Epetra_Map eInnerProdMap = toEpetra(InnerProdMap);
      RCP<Epetra_Vector> InnerProd_local = rcp(new Epetra_Vector(eInnerProdMap,true));

      if(left->getNodeNumRows() != right->getNodeNumRows()) std::cout << "Error: getNodeNumRows of left != right" << std::endl;

      for(size_t n=0; n<left->getNodeNumRows(); n++)
      {
        Teuchos::ArrayView<const LocalOrdinal> lindices_left;
        Teuchos::ArrayView<const Scalar> lvals_left;
        left->getLocalRowView(n, lindices_left, lvals_left);

        Teuchos::ArrayView<const LocalOrdinal> lindices_right;
        Teuchos::ArrayView<const Scalar> lvals_right;
        right->getLocalRowView(n, lindices_right, lvals_right);

        std::vector<Scalar> retvals;
        std::vector<GlobalOrdinal> retindx;
        retvals.reserve(lindices_left.size());  // we have not more than nnz_left possible matchings of left and right gids
        retindx.reserve(lindices_left.size());

        for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
        {
          for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
          {
            GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
            GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
            if(left_gid == right_gid)
            {
              retvals.push_back(lvals_left[i]*lvals_right[j]);
              retindx.push_back(left_gid); // note: we save the gids
              break; // skip remaining gids of right operator
            }
          }
        }

        InnerProd_local->SumIntoGlobalValues(Teuchos::as<size_t>(retvals.size()),&retvals[0],&retindx[0]); // we use the gids
      }

      // safety check
      if(InnerProd_local->GlobalLength() != InnerProd_local->MyLength()) std::cout << "Global and Local length have to be the same!" << std::endl;

      /////////////// extract all local gids and values
      std::vector<Scalar>           localvalues(InnerProd_local->GlobalLength(),Teuchos::ScalarTraits<Scalar>::zero());
      std::vector<GlobalOrdinal>    globalindices(InnerProd_local->GlobalLength(),-Teuchos::ScalarTraits<Scalar>::one());
      std::vector<Scalar>           globalvalues(InnerProd_local->GlobalLength(),Teuchos::ScalarTraits<Scalar>::zero());
      for(size_t i=0; i<localvalues.size(); i++)
      {
        localvalues[i] = (*InnerProd_local)[i];
        globalindices[i] = InnerProd_local->Map().GID(i); //i; // store indices (attention: local index = global index)
      }

      /////////////// sum up all values to global
      eInnerProdMap.Comm().SumAll(&localvalues[0],&globalvalues[0],localvalues.size()); // sum up all local vectors

      /////////////// save all global information in local replicated vector InnerProd_local
      InnerProd_local->ReplaceGlobalValues(globalvalues.size(),&globalvalues[0],&globalindices[0]);

      /////////////// free some memory
      localvalues.clear();
      globalindices.clear();
      globalvalues.clear();

      RCP<Xpetra::EpetraVector> xInnerProd_local = Teuchos::rcp(new Xpetra::EpetraVector(InnerProd_local));
      return xInnerProd_local;
      #else
      throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }

    // Tpetra-specific code
    else if(left->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      // Tpetra specific!
      RCP<Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > InnerProd_local =
          //Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(InnerProdMap, true);
          rcp(new Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(InnerProdMap,true));

      for(size_t n=0; n<left->getNodeNumRows(); n++)
      {
        Teuchos::ArrayView<const LocalOrdinal> lindices_left;
        Teuchos::ArrayView<const Scalar> lvals_left;

        left->getLocalRowView(n, lindices_left, lvals_left);


        Teuchos::ArrayView<const LocalOrdinal> lindices_right;
        Teuchos::ArrayView<const Scalar> lvals_right;

        right->getLocalRowView(n, lindices_right, lvals_right);


        for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
        {
          for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
          {
            GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
            GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
            if(left_gid == right_gid)
            {
              InnerProd_local->sumIntoLocalValue(left_gid, 0, lvals_left[i] * lvals_right[j]);
            }
          }
        }
      }

      InnerProd_local->reduce();

      return InnerProd_local;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }

    return Teuchos::null;
  }

  Teuchos::RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > BuildLocalReplicatedColMap(const RCP<Operator>& DinvAP0) const
  {
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
    Teuchos::RCP< const Teuchos::Comm< int > > comm = DinvAP0->getRangeMap()->getComm();

    // generate allreduced column-based array vector for local damping factors
    Teuchos::Array<GlobalOrdinal> localreplicatedcolgids;
    reduceAllXpetraMap(localreplicatedcolgids,*(DinvAP0->getDomainMap()));

    // Epetra-specific code
    if (DinvAP0->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      Teuchos::RCP<Epetra_Map> epetra_locreplicatedomegagids = Teuchos::rcp(new Epetra_Map(localreplicatedcolgids.size(),localreplicatedcolgids.size(),&localreplicatedcolgids[0],0,*Xpetra::toEpetra(comm)));
      const RCP< const Xpetra::EpetraMap> xpetra_locreplicatedomegagids = Teuchos::rcp(new Xpetra::EpetraMap(epetra_locreplicatedomegagids));
      return xpetra_locreplicatedomegagids;
#else
      throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }
    // Tpetra-specific code
    else if(DinvAP0->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      Teuchos::RCP<Teuchos::ArrayView<const GlobalOrdinal> > arView = Teuchos::rcp(new Teuchos::ArrayView<const GlobalOrdinal>(&localreplicatedcolgids[0],localreplicatedcolgids.size()));
      //Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tpetra_locreplicatedomegagids = Tpetra::createNonContigMap<GlobalOrdinal,LocalOrdinal>(*arView, comm);
      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tpetra_locreplicatedomegagids = Tpetra::createLocalMap<GlobalOrdinal,LocalOrdinal>(Teuchos::as<size_t>(DinvAP0->getColMap()->getMaxAllGlobalIndex()+1), comm); // create local map
      const RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > xpetra_locreplicatedomegagids = Xpetra::toXpetra(tpetra_locreplicatedomegagids);
      xpetra_locreplicatedomegagids->describe(*fos,Teuchos::VERB_EXTREME);
      return xpetra_locreplicatedomegagids;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }

    return Teuchos::null;
  }

  void Test(RCP<Operator> AP0) const {
    Teuchos::RCP< const Teuchos::Comm< int > > comm = AP0->getRangeMap()->getComm();
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

    if (AP0->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      RCP<const Epetra_CrsMatrix> epA = MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2EpetraCrs(AP0);
      std::cout << "todo" << std::endl;
#else
      throw(Exceptions::RuntimeError("Error."));
#endif
    } else if(AP0->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(AP0);


      RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > localomegamap =
          Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal, Node>::createLocalMap(Xpetra::UseTpetra, AP0->getColMap()->getMaxAllGlobalIndex(), comm);

      RCP<Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > localomegaVec =
          //Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(localomegamap);
          rcp(new Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(localomegamap,true));

      RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > omegaImport =
          Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(AP0->getColMap(),localomegamap);

      RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > distrVec =
          Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(AP0->getColMap());

      distrVec->putScalar(comm->getRank()+13);
      localomegaVec->putScalar(0.0);


      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;

      //AP0->describe(*fos,Teuchos::VERB_EXTREME);

      LocalOrdinal locOrd = AP0->getRowMap()->getLocalElement(52);
      std::cout << "PROC: " << comm->getRank() << " globalElement: " << 52 << " localElement: " << locOrd << std::endl;

      if(locOrd!=Teuchos::OrdinalTraits<LocalOrdinal>::invalid())
      {
        AP0->getLocalRowView(locOrd, lindices_right, lvals_right);
        for (LocalOrdinal i=0; i<lindices_right.size(); i++)
        {
          localomegaVec->replaceLocalValue(AP0->getColMap()->getGlobalElement(lindices_right[i]),0,lvals_right[i]);
          std::cout << lindices_right << " GIDS: " << AP0->getColMap()->getGlobalElement(lindices_right[i]) << std::endl;
          std::cout << lvals_right << std::endl;
        }
      }

      //localomegaVec->doImport(*distrVec,*omegaImport,Xpetra::ADD);
      localomegaVec->reduce();

      //localomegaVec->describe(*fos,Teuchos::VERB_EXTREME);

      //distrVec->describe(*fos,Teuchos::VERB_EXTREME);

      std::vector<GlobalOrdinal> gids;
      gids.push_back(25*(comm->getRank()+1));
      gids.push_back(30*(comm->getRank()+1));
      gids.push_back(35*(comm->getRank()+1));
      gids.push_back(40*(comm->getRank()+1));
      RCP<Teuchos::ArrayView<const GlobalOrdinal> > arView = rcp(new Teuchos::ArrayView<const GlobalOrdinal>(&gids[0],gids.size()));
      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tmap = Tpetra::createNonContigMap<GlobalOrdinal,LocalOrdinal>(*arView, comm);
      tmap->describe(*fos,Teuchos::VERB_EXTREME);

#else
  throw(Exceptions::RuntimeError("Error."));
#endif
    }
  }
private:

  //! @name helper function for allreducing a Xpetra::Map
  //@{

  int FindMyPos(size_t nummyelements, const Teuchos::Comm<int>& comm) const
  {
    std::cout << "begin FindMyPos" << std::endl;
    const int myrank = comm.getRank();
    const int numproc = comm.getSize();

    std::vector<size_t> snum (numproc,0); // local vector of no of elements
    std::vector<size_t> rnum (numproc);   // gobal vector of no of elements on proc
    snum[myrank] = nummyelements;

    std::cout << "FindMyPose::reduce all" << std::endl;
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, numproc, &snum[0], &rnum[0]);
    std::cout << "FindMyPos:: after reduce all" << std::endl;

    return std::accumulate(&rnum[0], &rnum[myrank], 0);
  }

  void reduceAllXpetraMap(Teuchos::Array<GlobalOrdinal>& rredundant, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const
  {
    const int mynodepos = FindMyPos(map.getNodeNumElements(), *(map.getComm()));

    std::vector<GlobalOrdinal> sredundant(map.getGlobalNumElements(),0);

    Teuchos::ArrayView< const GlobalOrdinal > gids = map.getNodeElementList();
    std::copy(gids.getRawPtr(), gids.getRawPtr() + map.getNodeNumElements(), &sredundant[mynodepos]); // check me

    rredundant.resize(map.getGlobalNumElements());
    GlobalOrdinal numEle = map.getGlobalNumElements();
    std::cout << "reduce All in XpetraMap" << std::endl;
    Teuchos::reduceAll(*(map.getComm()),Teuchos::REDUCE_SUM, numEle ,&sredundant[0], &rredundant[0]);
    std::cout << "after reduce All in XpetraMap" << std::endl;
  }

  //@}


private:

  //! Input factories
  RCP<PFactory> initialPFact_;        //! Ptentative Factory
  RCP<SingleLevelFactoryBase> AFact_; //! A Factory

  //! Factory parameters
  std::string diagonalView_;
};

//! Friend print function.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, PgPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
  os << "Printing an PgPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_PGPFACTORY_SHORT

#endif /* MUELU_PGPFACTORY_HPP_ */
