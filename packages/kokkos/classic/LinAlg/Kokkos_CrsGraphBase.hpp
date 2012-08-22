//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_CRSGRAPHBASE_HPP
#define KOKKOS_CRSGRAPHBASE_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <numeric>

namespace Kokkos {

  /*! @class CrsGraphBase
      @brief An abstract base class providing a template for Kokkos-level sparse graph objects.
      @ingroup kokkos_crs_ops
     
      The Tpetra classes do not utilize this base class for interacting with Kokkos-level objects, and 
      the local sparse graph objects are there not required to inherit from this interface. However, 
      this class illustrates the methods that are required by Tpetra objects, and therefore provides a 
      potential starting point for extending Tpetra via new local graph and matrix types.
      
      \tparam Ordinal Defines the type of the column indices;
                      same as the LocalOrdinal template parameter of the encapsulating Tpetra::CrsGraph object.
      \tparam Node    Kokkos Node type; same as the Node template parameter of the encapsulating Tpetra::CrsGraph object.
  */
  template <class Ordinal,
            class Node>
  class CrsGraphBase {
  public:
    typedef Ordinal  ordinal_type;
    typedef Node     node_type;

    //! @name Constructors/Destructor
    //@{

    //! Default constuctor.
    CrsGraphBase(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params);

    //! Destructor.
    virtual ~CrsGraphBase();

    //@}

    //! @name Data entry and accessor methods.
    //@{

    RCP<Node> getNode() const;

    //! Return the number of rows in the graph.
    Ordinal getNumRows() const;

    //! Return the number of columns in the graph.
    Ordinal getNumCols() const;

    /// \brief Submit the indices and offset for the graph.
    /**
          \pre indices for row \c r are inds[r], where \f$r \in [b,e)\f$, where \f$b = ptrs[r]\f$ and \f$e = ptrs[r-1]\f$
          \pre ptrs has getNumRows()+1 entries
          \pre ptrs[0] == 0
          \pre ptrs[getNumRows()] == inds.size()
     */
    virtual void setStructure(const ArrayRCP<const size_t> &ptrs,
                              const ArrayRCP<const Ordinal> &inds) = 0;

    //@}

  private:
    RCP<Node> node_;
    Ordinal numRows_, numCols_;
    RCP<ParameterList> params_;
  };


  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraphBase<Ordinal,Node>::CrsGraphBase(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params) 
  : node_(node)
  , numRows_(numRows)
  , numCols_(numCols)
  , params_(params)
  {}

  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraphBase<Ordinal,Node>::~CrsGraphBase() {
  }

  // ======= numrows ===========
  template <class Ordinal, class Node>
  Ordinal CrsGraphBase<Ordinal,Node>::getNumRows() const {
    return numRows_;
  }

  // ======= numcols ===========
  template <class Ordinal, class Node>
  Ordinal CrsGraphBase<Ordinal,Node>::getNumCols() const {
    return numCols_;
  }

  // ======= node ===========
  template <class Ordinal, class Node>
  RCP<Node> CrsGraphBase<Ordinal,Node>::getNode() const {
    return node_;
  }

} // namespace Kokkos

#endif /* KOKKOS_CRSGRAPHBASE_HPP */
