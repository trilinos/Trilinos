// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Tsqr_MessengerBase_hpp
#define __TSQR_Tsqr_MessengerBase_hpp

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  /// \class MessengerBase
  ///
  /// Interface for an object that performs collective communication.
  /// Each message contains some number of objects of scalar type
  /// Datum.  Datum must have a default constructor and a copy
  /// constructor, and taking its address must make sense (in terms of
  /// extracting the useful data).
  template< class Datum >
  class MessengerBase {
  public:
    virtual void 
    send (const Datum sendData[], 
	  const int sendCount, 
	  const int destProc, 
	  const int tag) = 0;

    virtual void 
    recv (Datum recvData[], 
	  const int recvCount, 
	  const int srcProc, 
	  const int tag) = 0;

    virtual void 
    swapData (const Datum sendData[], 
	      Datum recvData[], 
	      const int sendRecvCount, 
	      const int destProc, 
	      const int tag) = 0;
    
    virtual Datum 
    globalSum (const Datum& inDatum) = 0;

    virtual void
    globalVectorSum (const Datum inData[], 
		     Datum outData[], 
		     const int count) = 0;

    ///
    /// Assumes that Datum objects are less-than comparable by the
    /// underlying communication protocol.
    ///
    virtual Datum 
    globalMin (const Datum& inDatum) = 0;

    ///
    /// Assumes that Datum objects are less-than comparable by the
    /// underlying communication protocol.
    ///
    virtual Datum 
    globalMax (const Datum& inDatum) = 0;

    virtual void
    broadcast (Datum data[], 
	       const int count,
	       const int root) = 0;

    virtual int rank () const = 0;
    virtual int size () const = 0;
    virtual void barrier () const = 0;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_MessengerBase_hpp

