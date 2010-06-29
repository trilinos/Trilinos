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

    virtual int rank () const = 0;
    virtual int size () const = 0;
    virtual void barrier () const = 0;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_MessengerBase_hpp

