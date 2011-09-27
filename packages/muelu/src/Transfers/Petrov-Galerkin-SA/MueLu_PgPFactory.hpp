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
        coarseLevel.Request("P",initialPFact_.get());
        fineLevel.Request("A",AFact_.get());
    };
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
      */
    bool Build(Level& fineLevel, Level &coarseLevel) const {

        std::ostringstream buf; buf << coarseLevel.GetLevelID();
        RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP_"+buf.str()));
        timer->start(true);

        // Level Get
        RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_.get());
        RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_.get());

        if(restrictionMode_)
            A = Utils2<Scalar,LocalOrdinal,GlobalOrdinal>::Transpose(A,true); // build transpose of A explicitely

        fineLevel.Release("A", AFact_.get());
        coarseLevel.Release("P", initialPFact_.get());

        std::cout << "after Ptent has been calculated" << std::endl;

        // calculate D^{-1} * A * Ptent

        bool doFillComplete=true;
        bool optimizeStorage=false;
        RCP<Operator> DinvAP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

        doFillComplete=true;
        optimizeStorage=false;
        Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(A);
        Utils::MyOldScaleMatrix(DinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

        // calculate local damping factors

        // compute D^{-1} * A * D^{-1} * A * P0
        doFillComplete=true;
        optimizeStorage=false;
        RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAP0,false,doFillComplete,optimizeStorage);

        doFillComplete=true;
        optimizeStorage=false;
        Utils::MyOldScaleMatrix(DinvADinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

        // generate allreduced column-based array vector for local damping factors
        Teuchos::Array<GlobalOrdinal> localreplicatedcolgids;
        reduceAllXpetraMap(localreplicatedcolgids,*(DinvAP0->getDomainMap()));

        std::cout << localreplicatedcolgids.toString() << std::endl;

        Teuchos::RCP< const Teuchos::Comm< int > > comm = DinvADinvAP0->getRangeMap()->getComm();
        RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > map =
            Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal, Node>::createLocalMap(Xpetra::UseTpetra, 20, comm);

        Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
        //map->describe(*fos,Teuchos::VERB_EXTREME);

        RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > locVec =
                Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(map);
        //locVec->describe(*fos,Teuchos::VERB_EXTREME);


///////////////////

        RCP<const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rcolVec =
                Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvADinvAP0->getColMap());
        RCP<const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rdomVec =
                Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvADinvAP0->getDomainMap());

        RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > coldomImport =
                Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(DinvADinvAP0->getDomainMap(),DinvADinvAP0->getColMap());

        Teuchos::ArrayView< const GlobalOrdinal > row_globlist = DinvADinvAP0->getRowMap()->getNodeElementList();
        std::vector<LocalOrdinal> row_nodvec(row_globlist.size(),0);
        std::vector<LocalOrdinal> row_locvec(row_globlist.size(),0);
        Teuchos::ArrayView< LocalOrdinal > row_nodlist(row_nodvec);
        Teuchos::ArrayView< LocalOrdinal > row_loclist(row_locvec);
        DinvADinvAP0->getRowMap()->getRemoteIndexList(row_globlist,row_nodlist,row_loclist);
        std::cout << "rowgloblist " <<  row_globlist << std::endl;
        std::cout << "rownodlist " <<  row_nodlist << std::endl;
        std::cout << "rowloclist " <<  row_loclist << std::endl;

        Teuchos::ArrayView< const GlobalOrdinal > globlist = DinvADinvAP0->getColMap()->getNodeElementList();
        std::vector<LocalOrdinal> nodvec(globlist.size(),0);
        std::vector<LocalOrdinal> locvec(globlist.size(),0);
        Teuchos::ArrayView< LocalOrdinal > nodlist(nodvec);
        Teuchos::ArrayView< LocalOrdinal > loclist(locvec);
        DinvADinvAP0->getColMap()->getRemoteIndexList(globlist,nodlist,loclist);
        std::cout << "colgloblist " <<  globlist << std::endl;
        std::cout << "colnodlist " <<  nodlist << std::endl;
        std::cout << "colloclist " <<  loclist << std::endl;

////////////////////

        Teuchos::ArrayView<const LocalOrdinal> lindices_right;
        Teuchos::ArrayView<const Scalar> lvals_right;

        DinvADinvAP0->getLocalRowView(49, lindices_right, lvals_right);
        std::cout << "localindices for row l 49 on proc " << comm->getRank() << " " <<  lindices_right << std::endl;
        for(size_t i = 0; i< lindices_right.size(); i++)
            std::cout << "PROC " << comm->getRank() << " LID: " << loclist[lindices_right[i]] << " GID: " << DinvADinvAP0->getColMap()->getGlobalElement(loclist[lindices_right[i]]) << std::endl;
        std::cout << "localvals for row l 49 on proc " << comm->getRank() << " " <<  lvals_right << std::endl;

        ///////////////// minimize with respect to the (D^{-1} A)' D^{-1} A norm.
        //
        //               diag( R0 (A D^{-1}' D^{-1} A' D^{-1} A' R0' )
        //  omega = ---------------------------------------------------------
        //           diag( R0 A D^{-1} A D^{-1} D^{-1} A' D^{-1} A' R0' )
        //Teuchos::Array<Scalar> Numerator(localreplicatedcolgids.size(),0.0);
        //MultiplyAll(DinvAP0,DinvADinvAP0,Numerator);  // -> DinvAP0_subset

        //std::cout << Numerator << std::endl;
        // smooth P

        // do wiggle-cols handling?

        return true;
    }

    bool BuildP(Level &fineLevel, Level &coarseLevel) const {
        std::cout << "TODO: remove me" << std::endl;
    }

    //@}

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

    void MultiplyAll(const RCP<Operator>& left, const RCP<Operator>& right, Teuchos::Array<Scalar>& InnerProd) const
    {
        if(!left->getDomainMap()->isSameAs(*(right->getDomainMap()))) std::cout << "domain maps of left and right do not match" << std::endl;
        if(!left->getColMap()->isSameAs(*(right->getColMap()))) std::cout << "col maps of left and right do not match" << std::endl;
        if(!left->getRowMap()->isSameAs(*(right->getRowMap()))) std::cout << "row maps of left and right do not match" << std::endl;

        Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
        left->getDomainMap()->describe(*fos,Teuchos::VERB_EXTREME);

        // communicate InnerProd local
        Teuchos::RCP< const Teuchos::Comm< int > > comm = left->getRangeMap()->getComm();



#if 0
        Teuchos::Array<Scalar> InnerProd_local(InnerProd.size(),0.0);

        for (LocalOrdinal n=0; n<left->getNodeNumRows(); n++)
        {
            std::cout << "begin loop" << std::endl;
            GlobalOrdinal grid_left = left->getRowMap()->getGlobalElement(n);
            std::cout << "PROC: " << comm->getRank() << " global grid (left): " << grid_left << " with " << left->getNumEntriesInLocalRow(n) << " entries" << std::endl;

            //size_t nnz_left = left->getNumEntriesInLocalRow(n);

            Teuchos::ArrayView<const LocalOrdinal> gindices_left;
            Teuchos::ArrayView<const Scalar> gvals_left;

            left->getLocalRowView(n, gindices_left, gvals_left);
            std::cout << "local indices (left): " << gindices_left << std::endl;

            GlobalOrdinal grid_right = right->getRangeMap()->getGlobalElement(n);
            std::cout << "PROC: " << comm->getRank() << " global grid (right): " << grid_right << std::endl;

            //size_t nnz_right = right->getNumEntriesInLocalRow(n);

            Teuchos::ArrayView<const LocalOrdinal> gindices_right;
            Teuchos::ArrayView<const Scalar> gvals_right;

            right->getLocalRowView(n, gindices_right, gvals_right);
            std::cout << "local indices (right): " << gindices_right << std::endl;

            //if(gindices_right != gindices_left) std::cout << "oops. gindices_right != gindices_left" << std::endl;


            for(size_t i=0; i<gvals_left.size(); i++)
            {
                for(size_t j=0; j<gvals_right.size(); j++)
                {
                    if (gindices_left[i] == gindices_right[j])
                    {

                        //GlobalOrdinal ind = left->getDomainMap()->getGlobalElement(gindices_left[i]);
                        //std::cout << "globalcolindex (left): " << ind << std::endl;
                        InnerProd_local[gindices_left[i]] += gvals_left[i]*gvals_right[j];
                        break;
                    }
                }
            }

        }



        Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<GlobalOrdinal>(InnerProd_local.size()), &InnerProd_local[0], &InnerProd[0]);

#endif
    }

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
