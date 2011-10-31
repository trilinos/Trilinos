#ifndef MUELU_GENERICRFACTORY_DECL_HPP
#define MUELU_GENERICRFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

/*
 * MueLu_GenericRFactory.hpp
 *
 *  Created on: 20.09.2011
 *      Author: tobias
 */

#include <iostream>

#include "Xpetra_CrsOperator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class GenericRFactory class.
    @brief Factory for building restriction operators using a prolongator factory
  */

  template < class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps >
  class GenericRFactory : public RFactory {

#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    GenericRFactory(RCP<PFactory> PFact = Teuchos::null)
    ;

    //! Destructor.
    virtual ~GenericRFactory() ;
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const ;

    //@}

    //! @name Build methods.
    //@{

    void Build(Level & fineLevel, Level & coarseLevel) const ;

    void BuildR(Level & fineLevel, Level & coarseLevel) const ; //BuildR

    //@}

  private:

    //! P Factory
    RCP<PFactory> PFact_;

  }; //class TransPFactory

} //namespace MueLu

#define MUELU_GENERICRFACTORY_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_GENERICRFACTORY_DECL_HPP
