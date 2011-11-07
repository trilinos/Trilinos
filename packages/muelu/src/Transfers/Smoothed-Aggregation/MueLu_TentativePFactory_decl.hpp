#ifndef MUELU_TENTATIVEPFACTORY_DECL_HPP
#define MUELU_TENTATIVEPFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_LAPACK.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_TwoLevelFactoryBase.hpp" // TODO: inheritence of TentativePFactory
#include "MueLu_Level.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class TentativePFactory class.
    @brief Factory for building tentative prolongator.

    Factory for creating tentative prolongator.   Nullspace vectors are split across aggregates so that they
    have local support.  The vectors with local support are factored via LAPACK QR.  The Q becomes the
    tentative prolongator, and the R becomes the coarse nullspace. 
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class TentativePFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{
    
    /*! @brief Constructor.
      \param AggregationFact -- (optional) factory that creates aggregates.
    */
    TentativePFactory(RCP<const FactoryBase> aggregatesFact = Teuchos::null, RCP<const FactoryBase> nullspaceFact = Teuchos::null, RCP<const FactoryBase> AFact = Teuchos::null);
    
    //! Destructor.
    virtual ~TentativePFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level & fineLevel, Level & coarseLevel) const;

    //@}

    //! @name Set/Get methods.
    //@{

    void TentativeWithQR(bool value);

    bool TentativeWithQR();

    //@}

    //! @name Build methods.
    //@{

    void Build(Level & fineLevel, Level & coarseLevel) const;

    void BuildP(Level & fineLevel, Level & coarseLevel) const;

    //@}

  private:
    //! @name Static methods.
    //@{
    
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;
    
    /*! @brief Make tentative prolongator with QR.

      We note that the implementation would have been *much* easier
      if Ptent were allowed to have a row map based upon aggregates, i.e., a row
      map such that all DoF's in an aggregate were local and consecutive.
      In this case, we could still have used the matrix A's domain map to FillComplete Ptent.
      However, the prolongator smoothing step Ptent - A*Ptent would then be impossible
      because neither Epetra nor Tpetra support adding matrices with differing row maps.

      The following is a high-level view of how this method is implemented.  The result is a tentative
      prolongator such that Ptent has the same row map as A.

      1) The nullspace NS is communicated in such a way that if proc A owns aggregate k, the part of NS
         corresponding to k is replicated on A.  This will happen if aggregate k contains DoFs belonging
         to proc A and other procs.

      2) Each processor A does a local QR on the parts of the NS corresponding to aggregates that A owns.

      3) The rows of Q that A owns (i.e., rows of Q corresponding to DoFs that A owns) are inserted immediately into Ptentative.

      4) Any rows of Q that A doesn't own  are set  to the owning processor as follows:
         (This step is necessary because Epetra does not allow insertion of rows owned by other processors.)

         a) Q is stored as a CSR matrix.  We know each row has at most dim(NS) nonzeros.
            For each row of Q, we must transmit both the values and columns.  We do this by communicating separate vectors,
            each of length dim(NS), for the values and column indices.

         b) We try to make this communication more efficient by first communicating just the global row numbers
            that each processor should expect to receive. A "reduced" map is created from just the expected row numbers.

         c) We then communicate the rows of Q themselves using this "reduced" map, i.e., the target of the Importer is the reduced map.
            Otherwise, the only alternative was to base the Importer on A's rowmap, which is probably much larger than the reduced map.

      5) Once received, the rows are inserted by the owning processes and Ptent is fillCompleted.

      - FIXME There is no attempt to detect if the aggregate is too small to support the NS.
    */
     void MakeTentative(const Operator& fineA, const Aggregates& aggregates, const MultiVector & fineNullspace, //-> INPUT
                        RCP<MultiVector> & coarseNullspace, RCP<Operator> & Ptentative) const;                  //-> OUTPUT

    //@}

  private:
    RCP<const FactoryBase> aggregatesFact_; //! Factory that creates aggregates
    RCP<const FactoryBase> nullspaceFact_;  //! Factory creating the nullspace
    RCP<const FactoryBase> AFact_;          //! Define which matrix A is used in this factory

    bool QR_; //! use QR decomposition for improving nullspace information per default

  }; //class TentativePFactory

} //namespace MueLu

//TODO: noQR_

#define MUELU_TENTATIVEPFACTORY_SHORT
#endif // MUELU_TENTATIVEPFACTORY_DECL_HPP
