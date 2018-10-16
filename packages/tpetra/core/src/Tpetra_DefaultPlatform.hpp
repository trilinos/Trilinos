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

#ifndef TPETRA_DEFAULT_PLATFORM_HPP
#define TPETRA_DEFAULT_PLATFORM_HPP

#include <Kokkos_DefaultNode.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_SerialPlatform.hpp"
#ifdef HAVE_MPI
#  include "Tpetra_MpiPlatform.hpp"
#endif

namespace Tpetra {

/** \brief Returns a default platform appropriate for the enviroment.

    \warning This class is DEPRECATED and will be REMOVED SOON.  Do
      not use <tt>*Platform</tt> classes any more.  To initialize
      Tpetra, include <tt>Tpetra_Core.hpp</tt> and use
      Tpetra::ScopeGuard, or Tpetra::initialize and Tpetra::finalize.
      To get Tpetra's default Comm instance, include
      <tt>Tpetra_Core.hpp</tt> and call
      <tt>Tpetra::getDefaultComm()</tt>.  For the default Node type,
      use <tt>Tpetra::Map<>::node_type</tt>.  Do not create Node
      instances yourself.  It is OK for Node instances to be null.
 */
class TPETRA_DEPRECATED DefaultPlatform {
public:
  /// \brief The default platform type specified at compile time.
  ///
  /// \warning This typedef is DEPRECATED and will be removed soon!
#ifdef HAVE_TPETRA_MPI
  typedef MpiPlatform< ::Tpetra::Details::DefaultTypes::node_type> DefaultPlatformType;
#else
  typedef SerialPlatform< ::Tpetra::Details::DefaultTypes::node_type> DefaultPlatformType;
#endif

  /// \brief Return a reference to the default platform singleton.
  ///
  /// \warning This method is DEPRECATED and will be removed soon!
  static DefaultPlatformType& getDefaultPlatform ();

private:
  //! The default platform singleton.
  static Teuchos::RCP<DefaultPlatformType> platform_;
};

} // namespace Tpetra

#endif // TPETRA_DEFAULT_PLATFORM_HPP

