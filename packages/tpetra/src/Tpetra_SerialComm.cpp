/*Paul
04-Aug-2002 Status: Templated for class T. All Epetra methods except Distributor.
03-Sept-2002 Took out Directory and ImageID methods. Templated for PacketType, OrdinalType.
*/

namespace Tpetra {

  // default constructor
  template<class PacketType, class OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm() 
    : Object("Tpetra::Comm[Serial]") { cout << "(SC) SerialComm constructor called." << endl;}
  
  // copy constructor
  template<class PacketType, class OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm(const SerialComm<PacketType, OrdinalType>& Comm) 
    : Object(Comm.label()) {}
  
  // destructor
  template<class PacketType, class OrdinalType>
  SerialComm<PacketType, OrdinalType>::~SerialComm() { cout << "(SC) SerialComm destructor called." << endl; }
  
  // gather all
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      allVals[i] = myVals[i];
  }
  
  // sum
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      globalSums[i] = partialSums[i];
  }
  
  // min/max
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      globalMaxs[i] = partialMaxs[i];
  }
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      globalMins[i] = partialMins[i];
  }
  
  // scansum
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      scanSums[i] = myVals[i];
  }
  
} // namespace Tpetra
