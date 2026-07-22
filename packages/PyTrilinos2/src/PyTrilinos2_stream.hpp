#ifndef TEUCHOS_STREAM_HPP
#define TEUCHOS_STREAM_HPP

#include <iostream>
#include <fstream>
#include "Teuchos_RCP.hpp"

namespace Teuchos {

  std::ofstream openOfstream(std::string filename){
    std::ofstream outdata;
    outdata.open(filename);
    return outdata;
  }

  void closeOfstream(std::ofstream &outdata) {
    outdata.close();
  }

  Teuchos::RCP<std::ostream> getCout() {
    Teuchos::RCP<std::ostream> out = Teuchos::rcp(&std::cout, false);
    return out;
  }
}

#endif
