// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_SERIALPLATFORM_HPP
#define XPETRA_SERIALPLATFORM_HPP

#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_Describable.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Xpetra {

  //! \brief A implementation of the Platform class for serial platforms.
  template<class Node=KokkosClassic::DefaultNode::DefaultNodeType>
  class SerialPlatform : public Teuchos::Describable {
  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Node NodeType;
    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    explicit SerialPlatform(const Teuchos::RCP<Node> &node);

    //! Destructor
    ~SerialPlatform();

    //@}

    //! @name Class Creation and Accessor Methods
    //@{

    //! Comm Instance
    const Teuchos::RCP< const Teuchos::SerialComm<int> > getComm() const;

    //! Get Get a node for parallel computation.
    const Teuchos::RCP<Node> getNode() const;

    //@}
  private:
    SerialPlatform(const SerialPlatform<Node> &platform);

  protected:
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::SerialComm<int> > comm_;
  };

  template<class Node>
  SerialPlatform<Node>::SerialPlatform(const Teuchos::RCP<Node>& /* node */) :
    comm_ (Teuchos::rcp(new Teuchos::SerialComm<int>()))
  {}

  template<class Node>
  SerialPlatform<Node>::~SerialPlatform() {  }

  template<class Node>
  const Teuchos::RCP< const Teuchos::SerialComm<int> >
  SerialPlatform<Node>::getComm() const {
    return comm_;
  }

  template<class Node>
  const Teuchos::RCP< Node >
  SerialPlatform<Node>::getNode() const {
    return Teuchos::null;
  }

} // namespace Xpetra

#endif // XPETRA_SERIALPLATFORM_HPP
