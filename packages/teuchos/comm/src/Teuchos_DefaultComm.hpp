// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_DEFAULT_COMM_HPP
#define TEUCHOS_DEFAULT_COMM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif

namespace Teuchos {

/// \class DefaultComm
/// \brief Return a default global communicator appropriate for the build.
///
/// Use this class to get a \c Teuchos::Comm instance representing the
/// default global communicator.  If Teuchos was built with MPI (i.e.,
/// if the HAVE_MPI macro is defined), then the default communicator
/// wraps MPI_COMM_WORLD.  Otherwise, it is a "serial" communicator
/// (containing one process, whose rank is zero).
///
/// \tparam OrdinalType The ordinal type for the \c Comm communicator
///   wrapper template class.  \c Comm uses OrdinalType to represent
///   things like array lengths and indices.
///
/// \note (mfh 19 Jul 2011, 22 Dec 2011) OrdinalType is called
/// OrdinalType and not Ordinal, because of a bug in Intel's C++
/// compiler (version 11.1).  This compiler confuses the Ordinal
/// template parameter of DefaultComm with the Teuchos::Ordinal
/// typedef.  The Ordinal template parameter should actually shadow
/// the typedef in Teuchos, and it does with GCC 4.5.1, but does not
/// with Intel's compiler.  This may be the case with other compilers
/// as well, but I haven't tested them yet.  If you use Ordinal as the
/// template parameter in this class, the following line of code will
/// result in a compile error.
/// \code
/// RCP<const Comm<int> pComm = DefaultComm<int>::getDefaultSerialComm (null); 
/// \endcode
template<typename OrdinalType>
class DefaultComm {
public:

  /// \brief Return the default global communicator.
  ///
  /// \warning When running with MPI, do not call this function until
  ///   after MPI has been initialized!  You can use \c
  ///   GlobalMPISesssion to initialize MPI without explicitly
  ///   depending on the MPI interface or the mpi.h header file.  (If
  ///   Trilinos was not built with MPI, \c GlobalMPISession will do
  ///   the right thing, so you can use it unconditionally.)
  static Teuchos::RCP<const Comm<OrdinalType> > getComm();

  /// \brief Return a serial comm if the input comm is null.
  ///
  /// If the input communicator \c comm is null, return the default
  /// serial communicator.  Otherwise, just return the input.
  static Teuchos::RCP<const Comm<OrdinalType> >
  getDefaultSerialComm( const Teuchos::RCP<const Comm<OrdinalType> > &comm );

private:

  /// \brief The default global communicator.  
  ///
  /// If Teuchos was built with MPI, this is a wrapper for
  /// MPI_COMM_WORLD.  Otherwise, this is a "serial" communicator
  /// (containing one process, whose rank is zero).
  static Teuchos::RCP<const Comm<OrdinalType> > comm_;

  //! A "serial" communicator (containing one process, whose rank is zero).
  static Teuchos::RCP<const Comm<OrdinalType> > defaultSerialComm_;
};


template<typename OrdinalType>
Teuchos::RCP<const Teuchos::Comm<OrdinalType> >
DefaultComm<OrdinalType>::getComm()
{
  if(!comm_.get()) {
#ifdef HAVE_MPI
    comm_ = rcp(new MpiComm<OrdinalType>(opaqueWrapper((MPI_Comm)MPI_COMM_WORLD)));
#else // HAVE_MPI    
    comm_ = rcp(new SerialComm<OrdinalType>());
#endif // HAVE_MPI    
  }
  return comm_;
}

template<typename OrdinalType>
Teuchos::RCP<const Teuchos::Comm<OrdinalType> >
DefaultComm<OrdinalType>::getDefaultSerialComm(
  const Teuchos::RCP<const Comm<OrdinalType> > &comm
  )
{
  if( comm.get() )
    return comm;
  else
    return defaultSerialComm_;
}

template<typename OrdinalType>
Teuchos::RCP<const Teuchos::Comm<OrdinalType> >
DefaultComm<OrdinalType>::comm_ = Teuchos::null;

template<typename OrdinalType>
Teuchos::RCP<const Teuchos::Comm<OrdinalType> >
DefaultComm<OrdinalType>::defaultSerialComm_
= Teuchos::rcp(new Teuchos::SerialComm<OrdinalType>());

} // namespace Teuchos

#endif // TEUCHOS_DEFAULT_COMM_HPP
