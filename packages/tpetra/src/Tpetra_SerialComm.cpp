/*Paul
04-Aug-2002 Status: Templated for class T. All Epetra methods except Distributor.
*/

#include "Tpetra_SerialDirectory.h"

namespace Tpetra {

  // default constructor
  template<class PacketType, class OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm() 
    : Object("Tpetra::Comm[Serial]") {}
  
  // copy constructor
  template<class PacketType, class OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm(const SerialComm<PacketType, OrdinalType>& Comm) 
    : Object(Comm.label()) {}
  
  // destructor
  template<class PacketType, class OrdinalType>
  SerialComm<PacketType, OrdinalType>::~SerialComm() {}
  
  // gather all
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType count) const {
    for(OrdinalType i=0; i<count; i++)
      allVals[i] = myVals[i];
  }
  
  // sum
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType count) const {
    for(OrdinalType i=0; i<count; i++)
      globalSums[i] = partialSums[i];
  }
  
  // min/max
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType count) const {
    for(OrdinalType i=0; i<count; i++)
      globalMaxs[i] = partialMaxs[i];
  }
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType count) const {
    for(OrdinalType i=0; i<count; i++)
      globalMins[i] = partialMins[i];
  }
  
  // scansum
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType count) const {
    for(OrdinalType i=0; i<count; i++)
      scanSums[i] = myVals[i];
  }
  
  // directory
  template<class PacketType, class OrdinalType>
  Directory<OrdinalType>* SerialComm<PacketType, OrdinalType>::createDirectory(const ElementSpace<OrdinalType>& ElementSpace) const {
    Directory<OrdinalType>* dir = dynamic_cast<Directory<OrdinalType>*>(new SerialDirectory<OrdinalType>(ElementSpace));
    return(dir);
  }
  
  // print
  template<class PacketType, class OrdinalType>
  void SerialComm<PacketType, OrdinalType>::print(ostream& os) const {
    os << "::Memory Image " << getMyImageID() << " of " << getNumImages() << " total images";
  }
  
} // namespace Tpetra
