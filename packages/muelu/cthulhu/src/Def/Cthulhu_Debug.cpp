#include <iostream>
#include <Teuchos_Hashtable.hpp>
#include <Teuchos_Array.hpp>

#include "Cthulhu_Debug.hpp"

// TODO Description of this file

// h{string} = integer;
typedef Teuchos::Hashtable<std::string, int> htable;

// h{string}{string} = integer;
typedef Teuchos::Hashtable<std::string, htable> hhtable;

// debugMeHashTable{fileName}{funcName} = number of occurrence
hhtable debugMeHashTable;

//TODO: description
void cthulhu_debug_me(const std::string & file, const std::string & funcName) {
  hhtable & h = debugMeHashTable;

  // fileName = basename(file);
  int pos = file.rfind("/");
  std::string fileName = file.substr(pos+1);

  if (h.containsKey(fileName)) {
    htable & subH = const_cast<htable &>(h.get(fileName)); // const_cast :-(

    int nb=0;    
    if (subH.containsKey(funcName)) {
      nb = subH.get(funcName);
      // std::cout << nb << std::endl;
    }

    subH.put(funcName, nb+1);
    

  } else {

    htable subH;
    subH.put(funcName,1);
    h.put(fileName, subH);

  }

}

// TODO: description
void cthulhu_debug_me_print() {
  Teuchos::Hashtable< std::string, Teuchos::Hashtable<std::string, int > > & h = debugMeHashTable;
  
  // std::cout << h.toString() << std::endl;
  
  Teuchos::Array< std::string > keys;
  Teuchos::Array< htable > values;
  
  h.arrayify(keys, values);
  
  for(int i=0; i<keys.length(); i++) {
    std::cout << keys[i] << ": " << std::endl;

    Teuchos::Array< std::string > subKeys;
    Teuchos::Array< int > subValues;
    
    values[i].arrayify(subKeys, subValues);

    for(int j=0; j<subKeys.length(); j++) {
      std::cout << "  - " << subKeys[j] << " (" << subValues[j] << ")" << std::endl;
    }

  }
  
}
