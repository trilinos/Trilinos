// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef TEUCHOS_SERIAL_COMM_HPP
#define TEUCHOS_SERIAL_COMM_HPP

#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"

namespace Teuchos {

/** \brief Concrete serial communicator subclass.
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal>
class SerialComm : public Comm<Ordinal> {
public:

  //! @name Constructors 
  //@{

  /** \brief . */
  SerialComm();

  //@}

  //! @name Overridden from Comm 
  //@{

  /** \brief . */
  virtual int getRank() const;
  /** \brief . */
  virtual int getSize() const;
  /** \brief . */
  virtual void barrier() const;
  /** \brief . */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const;
  /** \brief . */
  virtual void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const;
  /** \brief . */
  virtual void reduceAllAndScatter(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvCounts[], char myGlobalReducts[]
    ) const;
  /** \brief . */
	virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> ireceive(
    const ArrayView<char> &Buffer,
    const int sourceRank
    ) const;
  /** \brief . */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest> > &requests
    ) const;
  /** \brief . */
  virtual void wait(
    const Ptr<RCP<CommRequest> > &request
    ) const;

  //@}

  //! @name Overridden from Describable 
  //@{

  /** \brief . */
  std::string description() const;

  //@}
	
};

// ////////////////////////
// Implementations

// Constructors

template<typename Ordinal>
SerialComm<Ordinal>::SerialComm()
{}

// Overridden from Comm
  
template<typename Ordinal>
int SerialComm<Ordinal>::getRank() const
{
  return 0;
}
  
template<typename Ordinal>
int SerialComm<Ordinal>::getSize() const
{
  return 1;
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::barrier() const
{
  // Nothing to do
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::broadcast(
  const int /*rootRank*/, const Ordinal /*bytes*/, char []/*buffer*/
  ) const
{
  // Nothing to do
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::gatherAll(
  const Ordinal sendBytes, const char sendBuffer[]
  ,const Ordinal recvBytes, char recvBuffer[]
  ) const
{
  (void)sendBytes;  // to remove "unused parameter" warning
  (void)recvBytes;
  (void)sendBuffer;
  (void)recvBuffer;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!(sendBytes==recvBytes));
#endif
  std::copy(sendBuffer,sendBuffer+sendBytes,recvBuffer);
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::reduceAll(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
  ) const
{
  (void)reductOp;
  std::copy(sendBuffer,sendBuffer+bytes,globalReducts);
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::reduceAllAndScatter(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal sendBytes, const char sendBuffer[]
  ,const Ordinal recvCounts[], char myGlobalReducts[]
  ) const
{
  // Ignore unused arguments
  (void)reductOp;
  (void)sendBytes;
  (void)sendBuffer;
  (void)recvCounts;
  (void)myGlobalReducts;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( recvCounts==NULL || recvCounts[0] != sendBytes ); 
#endif
  std::copy(sendBuffer,sendBuffer+sendBytes,myGlobalReducts);
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::scan(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
  ) const
{
  (void)reductOp;
  std::copy(sendBuffer,sendBuffer+bytes,scanReducts);
}
  
template<typename Ordinal>
void SerialComm<Ordinal>::send(
  const Ordinal /*bytes*/, const char []/*sendBuffer*/, const int /*destRank*/
  ) const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::send(...): Error, you can not call send(...) when you"
    " only have one process!"
    );
}
  
template<typename Ordinal>
int SerialComm<Ordinal>::receive(
  const int /*sourceRank*/, const Ordinal /*bytes*/, char []/*recvBuffer*/
  ) const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::receive(...): Error, you can not call receive(...) when you"
    " only have one process!"
    );
  // The next line will never be reached, but a return is required on some platforms
  return 0; 
}

template<typename Ordinal>
void SerialComm<Ordinal>::readySend(
  const ArrayView<const char> &/*sendBuffer*/,
  const int /*destRank*/
  ) const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::readySend(...): Error, you can not call readySend(...) when you"
    " only have one process!"
    );
}



template<typename Ordinal>
RCP<CommRequest> SerialComm<Ordinal>::isend(
  const ArrayView<const char> &/*sendBuffer*/,
  const int /*destRank*/
  ) const
{
  TEST_FOR_EXCEPT(true);
  return null;
}


template<typename Ordinal>
RCP<CommRequest> SerialComm<Ordinal>::ireceive(
  const ArrayView<char> &/*Buffer*/,
  const int /*sourceRank*/
  ) const
{
  TEST_FOR_EXCEPT(true);
  return null;
}


template<typename Ordinal>
void SerialComm<Ordinal>::waitAll(
  const ArrayView<RCP<CommRequest> > &/*requests*/
  ) const
{
  TEST_FOR_EXCEPT(true);
}


template<typename Ordinal>
void SerialComm<Ordinal>::wait(
  const Ptr<RCP<CommRequest> > &/*request*/
  ) const
{
  TEST_FOR_EXCEPT(true);
}


// Overridden from Describable

template<typename Ordinal>
std::string SerialComm<Ordinal>::description() const
{
  std::ostringstream oss;
  oss << "Teuchos::SerialComm<"<<OrdinalTraits<Ordinal>::name()<<">";
  return oss.str();
}

} // namespace Teuchos

#endif // TEUCHOS_SERIAL_COMM_HPP
