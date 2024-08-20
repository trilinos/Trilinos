// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_SERIAL_COMM_HPP
#define TEUCHOS_SERIAL_COMM_HPP

#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"


namespace Teuchos {

/// \class SerialCommStatus
/// \brief Implementation of \c CommStatus for a serial communicator.
///
/// \tparam OrdinalType The same template parameter as \c Comm.  Only
///   use \c int here.  We only make this a template class for
///   compatibility with \c Comm.
template<class OrdinalType>
class SerialCommStatus : public CommStatus<OrdinalType> {
public:
  //! Default constructor.
  SerialCommStatus () {}

  //! The source rank that sent the message (must be zero).
  OrdinalType getSourceRank () { return 0; }

  //! The tag of the received message.
  OrdinalType getTag () { return 0; }
};


/** \brief Concrete serial communicator subclass.
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal>
class SerialComm : public Comm<Ordinal> {
public:
  /// \brief The current tag.
  ///
  /// \warning This method is ONLY for use by Teuchos developers.
  ///   Users should not depend on the interface of this method.
  ///   It may change or disappear at any time without warning.
  int getTag () const { return 0; }

  //! @name Constructors
  //@{

  /** \brief . */
  SerialComm();

  /** \brief Default copy constructor. */
  SerialComm(const SerialComm<Ordinal>& other);

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
  //! Gather values from all processes to the root process.
  virtual void
  gather (const Ordinal sendBytes, const char sendBuffer[],
          const Ordinal recvBytes, char recvBuffer[],
          const int root) const;
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
        virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual void
  send (const Ordinal bytes,
        const char sendBuffer[],
        const int destRank,
        const int tag) const;
  /** \brief . */
  virtual void ssend(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual void
  ssend (const Ordinal bytes,
         const char sendBuffer[],
         const int destRank,
         const int tag) const;
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
  virtual void
  readySend (const Ordinal bytes,
             const char sendBuffer[],
             const int destRank,
             const int tag) const;
  /** \brief . */
  virtual RCP<CommRequest<Ordinal> > isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  //! Variant of isend() that takes a tag.
  virtual RCP<CommRequest<Ordinal> >
  isend (const ArrayView<const char> &sendBuffer,
         const int destRank,
         const int tag) const;
  /** \brief . */
  virtual RCP<CommRequest<Ordinal> > ireceive(
    const ArrayView<char> &Buffer,
    const int sourceRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest<Ordinal> >
  ireceive (const ArrayView<char> &Buffer,
            const int sourceRank,
            const int tag) const;
  /** \brief . */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest<Ordinal> > > &requests
    ) const;
  /** \brief . */
  virtual void
  waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
           const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const;
  /** \brief . */
  virtual RCP<CommStatus<Ordinal> >
  wait (const Ptr<RCP<CommRequest<Ordinal> > >& request) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > duplicate() const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > split(const int color, const int key) const;
  /** brief . */
  virtual RCP< Comm<Ordinal> > createSubcommunicator(
    const ArrayView<const int> & ranks) const;

  //@}

  //! @name Overridden from Describable
  //@{

  /** \brief . */
  std::string description() const;

  //@}

};


/** \brief Nonmember constructor.
 *
 * \relates SerialComm
 */
template<typename Ordinal>
RCP<SerialComm<Ordinal> > createSerialComm()
{
  return Teuchos::rcp(new SerialComm<Ordinal>);
}


// ////////////////////////
// Implementations


// Constructors


template<typename Ordinal>
SerialComm<Ordinal>::SerialComm()
{}

template<typename Ordinal>
SerialComm<Ordinal>::SerialComm(const SerialComm<Ordinal>& /* other */)
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
  TEUCHOS_TEST_FOR_EXCEPT(!(sendBytes==recvBytes));
#endif
  std::copy(sendBuffer,sendBuffer+sendBytes,recvBuffer);
}


template<typename Ordinal>
void
SerialComm<Ordinal>::gather (const Ordinal sendBytes,
                             const char sendBuffer[],
                             const Ordinal recvBytes,
                             char recvBuffer[],
                             const int root) const
{
  (void) sendBytes;  // to remove "unused parameter" warning
  (void) recvBytes;
  (void) sendBuffer;
  (void) recvBuffer;
  (void) root;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!(sendBytes==recvBytes));
#endif
  std::copy (sendBuffer, sendBuffer + sendBytes, recvBuffer);
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::send(...): Error, you can not call send(...) when you"
    " only have one process!"
    );
}

template<typename Ordinal>
void SerialComm<Ordinal>::
send (const Ordinal /*bytes*/,
      const char []/*sendBuffer*/,
      const int /*destRank*/,
      const int /*tag*/) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::send(...): Error, you can not call send(...) when you"
    " only have one process!"
    );
}

template<typename Ordinal>
void SerialComm<Ordinal>::ssend(
  const Ordinal /*bytes*/, const char []/*sendBuffer*/, const int /*destRank*/
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::send(...): Error, you can not call ssend(...) when you"
    " only have one process!"
    );
}

template<typename Ordinal>
void
SerialComm<Ordinal>::ssend (const Ordinal /* bytes */,
                            const char* /* sendBuffer */,
                            const int /* destRank */,
                            const int /* tag */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::send(...): Error, you can not call ssend(...) when you"
    " only have one process!"
    );
}

template<typename Ordinal>
int SerialComm<Ordinal>::receive(
  const int /*sourceRank*/, const Ordinal /*bytes*/, char []/*recvBuffer*/
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::receive(...): Error, you can not call receive(...) when you"
    " only have one process!"
    );
}


template<typename Ordinal>
void SerialComm<Ordinal>::readySend(
  const ArrayView<const char> &/*sendBuffer*/,
  const int /*destRank*/
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::readySend(...): Error, you can not call readySend(...) when you"
    " only have one process!"
    );
}

template<typename Ordinal>
void
SerialComm<Ordinal>::readySend (const Ordinal bytes,
                                const char sendBuffer[],
                                const int destRank,
                                const int tag) const
{
  (void) bytes;
  (void) sendBuffer;
  (void) destRank;
  (void) tag;

  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"SerialComm<Ordinal>::readySend(...): Error, you can not call readySend(...) when you"
    " only have one process!"
    );
}

template<typename Ordinal>
RCP<CommRequest<Ordinal> > SerialComm<Ordinal>::isend(
  const ArrayView<const char> &/*sendBuffer*/,
  const int /*destRank*/
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "SerialComm<Ordinal>::isend: You cannot call isend when you only have one process." );
}


template<typename Ordinal>
RCP<CommRequest<Ordinal> >
SerialComm<Ordinal>::
isend (const ArrayView<const char> &/*sendBuffer*/,
       const int /*destRank*/,
       const int /*tag*/) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "SerialComm<Ordinal>::isend: You cannot call isend when you only have one process." );
}


template<typename Ordinal>
RCP<CommRequest<Ordinal> > SerialComm<Ordinal>::ireceive(
  const ArrayView<char> &/*Buffer*/,
  const int /*sourceRank*/
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "SerialComm<Ordinal>::ireceive: You cannot call isend when you only have one process." );
}


template<typename Ordinal>
RCP<CommRequest<Ordinal> >
SerialComm<Ordinal>::
ireceive (const ArrayView<char> &/*Buffer*/,
          const int /*sourceRank*/,
          const int /*tag*/) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "SerialComm<Ordinal>::ireceive: You cannot call isend when you only have one process." );
}


template<typename Ordinal>
void SerialComm<Ordinal>::waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests) const
{
  (void) requests;
  // There's nothing to wait on!
}


template<typename Ordinal>
void
SerialComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
         const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(statuses.size() < requests.size(),
    std::invalid_argument, "Teuchos::SerialComm::waitAll: There are not enough "
    "entries in the statuses array to hold all the results of the communication"
    " requests.  requests.size() = " << requests.size() << " > statuses.size() "
    "= " << statuses.size() << ".");

  for (typename ArrayView<RCP<CommRequest<Ordinal> > >::iterator it = requests.begin();
       it != requests.end(); ++it) {
    *it = null; // A postcondition of the Teuchos::Comm interface.
  }
}

template<typename Ordinal>
RCP<CommStatus<Ordinal> >
SerialComm<Ordinal>::wait (const Ptr<RCP<CommRequest<Ordinal> > > & request) const
{
  (void) request;
  TEUCHOS_TEST_FOR_EXCEPTION(request.getRawPtr() == NULL, std::invalid_argument,
    "Teuchos::SerialComm::wait: On input, the request pointer is null.");

  if (is_null (*request)) {
    return null; // Nothing to wait on...
  }
  *request = null;
  return rcp (new SerialCommStatus<Ordinal>);
}

template< typename Ordinal>
RCP< Comm<Ordinal> >
SerialComm<Ordinal>::duplicate() const
{
  return rcp(new SerialComm<Ordinal>(*this));
}

template<typename Ordinal>
RCP< Comm<Ordinal> >
SerialComm<Ordinal>::split(const int color, const int /*key*/) const
{
  if (color < 0) {
    return RCP< Comm<Ordinal> >();
  }
  // Simply return a copy of this communicator.
  return rcp(new SerialComm<Ordinal>(*this));
}

template<typename Ordinal>
RCP< Comm<Ordinal> >
SerialComm<Ordinal>::createSubcommunicator(const ArrayView<const int> &ranks) const
{
  if ((ranks.size()) == 1 && (ranks[0] == 0)) {
    return rcp(new SerialComm<Ordinal>(*this));
  } else {
    return RCP< Comm<Ordinal> >();
  }
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
