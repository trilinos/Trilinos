/*Paul
04-Aug-2002 Status: Templated for class T. All Epetra methods except Distributor.
03-Sept-2002 Took out Directory and ImageID methods. Templated for PacketType, OrdinalType.
12-Oct-2002 Added some consts (still some left).
12-Nov-2002 Changed remaining template<class...> to template<typename...>
*/

namespace Tpetra {

  // default constructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm() 
    : Object("Tpetra::Comm[Serial]") {}
  
  // copy constructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm(const SerialComm<PacketType, OrdinalType>& Comm) 
    : Object(Comm.label()) {}
  
  // destructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::~SerialComm() {}
  
  // gather all
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::gatherAll(PacketType* myVals, PacketType* allVals, const OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      allVals[i] = myVals[i];
  }
  
  // sum
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::sumAll(PacketType* partialSums, PacketType* globalSums, const OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      globalSums[i] = partialSums[i];
  }
  
  // min/max
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::maxAll(PacketType* partialMaxs, PacketType* globalMaxs, const OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      globalMaxs[i] = partialMaxs[i];
  }
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::minAll(PacketType* partialMins, PacketType* globalMins, const OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      globalMins[i] = partialMins[i];
  }
  
  // scansum
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::scanSum(PacketType* myVals, PacketType* scanSums, const OrdinalType count) const {
    for(OrdinalType i=0; i < count; i++)
      scanSums[i] = myVals[i];
  }

	// print

  
} // namespace Tpetra
