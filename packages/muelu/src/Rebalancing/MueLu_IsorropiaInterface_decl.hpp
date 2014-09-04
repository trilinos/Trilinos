/*
 * MueLu_IsorropiaInterface_decl.hpp
 *
 *  Created on: Jun 10, 2013
 *      Author: tobias
 */

#ifndef MUELU_ISORROPIAINTERFACE_DECL_HPP_
#define MUELU_ISORROPIAINTERFACE_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

//#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp> //TODO

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsGraph.hpp>
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraCrsGraph.hpp>
#endif

#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

  /*!
    @class IsorropiaInterface
    @brief Interface to Isorropia package.

  */

  //FIXME: this class should not be templated
  template <class LocalOrdinal = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType,
            class LocalMatOps = void>
  class IsorropiaInterface : public SingleLevelFactoryBase {

    typedef double Scalar; // FIXME
#undef MUELU_ISORROPIAINTERFACE_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    //! Constructor
    IsorropiaInterface() { }

    //! Destructor
    virtual ~IsorropiaInterface() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! @name Input
    //@{
    void DeclareInput(Level & level) const;
    //@}

    //! @name Build methods.
    //@{
    void Build(Level &level) const;

    //@}



  private:



  };  //class IsorropiaInterface

} //namespace MueLu

#define MUELU_ISORROPIAINTERFACE_SHORT
//#endif //if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)


#endif /* MUELU_ISORROPIAINTERFACE_DECL_HPP_ */
