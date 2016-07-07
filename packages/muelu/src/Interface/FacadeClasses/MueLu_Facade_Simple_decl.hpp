// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DECL_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DECL_HPP_

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>

#include "MueLu_FacadeClassBase_decl.hpp"

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class FacadeSimple : public FacadeClassBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors
    //@{

    /*! @brief Constructor that accepts a user-provided ParameterList.

        Constructor for parameter list interpreter which directly interprets Teuchos::ParameterLists

        @details The parameter list can be either in the easy parameter list format or in the factory driven parameter list format.

        @param[in] paramList (Teuchos::ParameterList): ParameterList containing the MueLu parameters
        @param[in] comm  (RCP<Teuchos::Comm<int> >): Optional RCP of a Teuchos communicator  (default: Teuchos::null)
        @param[in] factFact  (RCP<FactoryFactory>): Optional parameter allowing to define user-specific factory interpreters for user-specific extensions of the XML interface. (default: Teuchos::null)

     */
    FacadeSimple();

    //! Destructor.
    virtual ~FacadeSimple() { }

    //@}

    /*! @brief Set parameter list for FacadeClassFactory interpreter.

       The routine checks whether it is a parameter list in the easy parameter format or the more advanced factory-based parameter format and calls the corresponding interpreter routine.

       When finished, the parameter list is set that will used by the hierarchy build phase.

       This method includes validation and some pre-parsing of the list for:
           - verbosity level
           - data to export
           - cycle type
           - max coarse size
           - max levels
           - number of equations

       @param[in] paramList: ParameterList containing the MueLu parameters.
    */
    Teuchos::RCP<Teuchos::ParameterList> SetParameterList(const Teuchos::ParameterList& paramList);

    //! Call the SetupHierarchy routine from the HiearchyManager object.
    //void SetupHierarchy(Hierarchy& H) const;

  private:

    //! @brief String equivalent of preconditioner layout
    //static const std::string stringTemplate_;

    //! @ brief String defining default parameter list for facade input
    //static const std::string defaultParams_;
  };

} // namespace MueLu



#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DECL_HPP_ */
