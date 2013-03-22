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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_DEFAULTPLATFORM_HPP
#define XPETRA_DEFAULTPLATFORM_HPP

#include <Kokkos_DefaultNode.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_SerialPlatform.hpp"
#ifdef HAVE_MPI
#  include "Xpetra_MpiPlatform.hpp"
#endif

#include "Xpetra_ConfigDefs.hpp"

namespace Xpetra {

  /** \brief Returns a default platform appropriate for the enviroment.

  The DefaultPlatform mechanism is useful for easily accessing default
  Comm and Node types on a particular system.

  If HAVE_MPI is defined, then an instance of <tt>MpiPlatform</tt> will be
  created.  Otherwise, a <tt>SerialPlatform</tt> is returned.
  */
  class DefaultPlatform {
  public:
    //! Typedef indicating the default platform type specified at compile time. For a serial build, this will be SerialPlatform. Otherwise, it will be MpiPlatform.
#ifdef HAVE_MPI
    typedef MpiPlatform<Kokkos::DefaultNode::DefaultNodeType> DefaultPlatformType;
#else
    typedef SerialPlatform<Kokkos::DefaultNode::DefaultNodeType> DefaultPlatformType;
#endif

    /** \brief Return the default platform.
     */
    static DefaultPlatformType & getDefaultPlatform();

  private:

    static Teuchos::RCP<DefaultPlatformType> platform_;

  };

} // namespace Xpetra

#endif // XPETRA_DEFAULTPLATFORM_HPP
