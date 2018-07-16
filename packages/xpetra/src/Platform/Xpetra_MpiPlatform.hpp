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
#ifndef XPETRA_MPIPLATFORM_HPP
#define XPETRA_MPIPLATFORM_HPP

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"

namespace Xpetra {

  //! \brief A implementation of the Platform class for MPI-based platforms.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal.
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template<class Node=KokkosClassic::DefaultNode::DefaultNodeType>
  class MpiPlatform : public Teuchos::Describable {
  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Node NodeType;
    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    explicit MpiPlatform(Teuchos::RCP<Node> node);

    //! Constructor
    MpiPlatform(Teuchos::RCP<Node> node, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm);

    //! Destructor
    ~MpiPlatform();

    //@}

    //! @name Class Creation and Accessor Methods
    //@{

    //! Comm Instance
    Teuchos::RCP< const Teuchos::Comm<int> > getComm() const;

    //! Get a Kokkos Node instance.
    Teuchos::RCP<Node> getNode() const;

    //@}

  private:
    Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
    MpiPlatform(const MpiPlatform<Node> &platform);
  };

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> /* node */, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm) :
    comm_ (Teuchos::createMpiComm<int>(rawMpiComm))
  {}

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> /* node */) :
    comm_ (Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD)))
  {}

  template <class Node>
  MpiPlatform<Node>::~MpiPlatform() {  }

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(const MpiPlatform<Node> &platform) {
    comm_ = platform.comm_;
  }

  template <class Node>
  Teuchos::RCP< const Teuchos::Comm<int> >
  MpiPlatform<Node>::getComm() const {
    return comm_;
  }

  template <class Node>
  Teuchos::RCP<Node> MpiPlatform<Node>::getNode() const
  {  return Teuchos::null; }

} // namespace Xpetra

#endif // XPETRA_MPIPLATFORM_HPP

