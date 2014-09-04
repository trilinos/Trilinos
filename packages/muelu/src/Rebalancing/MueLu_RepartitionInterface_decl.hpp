/*
 * MueLu_RepartitionInterface_decl.hpp
 *
 *  Created on: 5 Sep 2013
 *      Author: wiesner
 */

#ifndef MUELU_REPARTITIONINTERFACE_DECL_HPP_
#define MUELU_REPARTITIONINTERFACE_DECL_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"


namespace MueLu {

  /*!
    @class RepartitionInterface
    @brief Helper class which transforms an "AmalgamatedPartition" array to an unamalgamated "Partition".

    This class is meant to be used with IsorropiaInterface which in general provides the amalgamated partition information
    and an AmalgamationFactory which defines the amalgamation/unamalgamation process.
    The output is a "Partition" (unamalgamated) which can be used by the RepartitionFactory class.

    Input: matrix A, unamalgamation information (that corresponds to matrix A)
  */

  //FIXME: this class should not be templated
  template <class LocalOrdinal = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class RepartitionInterface : public SingleLevelFactoryBase {

    typedef double Scalar; // FIXME
#undef MUELU_REPARTITIONINTERFACE_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    //! Constructor
    RepartitionInterface() { }

    //! Destructor
    virtual ~RepartitionInterface() { }
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



  };  //class RepartitionInterface

} //namespace MueLu

#define MUELU_REPARTITIONINTERFACE_SHORT
#endif /* MUELU_REPARTITIONINTERFACE_DECL_HPP_ */
