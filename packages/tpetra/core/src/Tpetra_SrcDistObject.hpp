// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_SRCDISTOBJECT_DECL_HPP
#define TPETRA_SRCDISTOBJECT_DECL_HPP

/// \file Tpetra_SrcDistObject.hpp
/// \brief Abstract base class for sources of an Import or Export.

#include <Teuchos_RCP.hpp>

namespace Tpetra {
  /// \class SrcDistObject
  /// \brief Abstract base class for objects that can be the source of
  ///   an Import or Export operation.
  ///
  /// Any object that may be the source of an Import or Export data
  /// redistribution operation must inherit from this class.  This
  /// class implements no methods, other than a trivial virtual
  /// destructor.  If a subclass X inherits from this class, that
  /// indicates that the subclass can be the source of an Import or
  /// Export, for <i>some set</i> of subclasses of DistObject.  A
  /// subclass Y of DistObject which is the target of the Import or
  /// Export operation will attempt to cast the input source
  /// SrcDistObject to a subclass which it knows how to treat as a
  /// source object.  The target subclass Y is responsible for knowing
  /// what source classes to expect, and how to interpret the
  /// resulting source object.
  ///
  /// DistObject inherits from this class, since a DistObject subclass
  /// may be either the source or the target of an Import or Export.
  /// A SrcDistObject subclass which does not inherit from DistObject
  /// need only be a valid source of an Import or Export; it need not
  /// be a valid target.
  ///
  /// This object compares to the Epetra class Epetra_SrcDistObject.
  /// Unlike in Epetra, this class in Tpetra does <tt>not</tt> include
  /// a getMap() method.  This is for two reasons.  First, consider
  /// the following inheritance hierarchy: DistObject and RowGraph
  /// inherit from SrcDistObject, and CrsGraph inherits from
  /// DistObject and RowGraph.  If SrcDistObject had a virtual getMap
  /// method, that would make resolution of the method ambiguous.
  /// Second, it is not necessary for SrcDistObject to have a getMap
  /// method, because a SrcDistObject alone does not suffice as the
  /// source of an Import or Export.  Any DistObject subclass must
  /// cast the SrcDistObject to a subclass which it knows how to treat
  /// as the source of an Import or Export.  Thus, it's not necessary
  /// for SrcDistObject to have a getMap method, since it needs to be
  /// cast anyway before use.  In general, I prefer to keep interfaces
  /// as simple as possible.
  class SrcDistObject {
  public:
    /// \brief Virtual destructor.
    ///
    /// It's necessary to provide this both for memory safety of
    /// derived classes, and to make SrcDistObject a polymorphic type
    /// (otherwise it can't be the argument of a dynamic_cast).
    virtual ~SrcDistObject () {};
  };

} // namespace Tpetra

#endif /* TPETRA_SRCDISTOBJECT_DECL_HPP */
