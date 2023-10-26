#ifndef PYROL_STREAM
#define PYROL_STREAM

#include <iostream>
#include <fstream>

namespace ROL {

  std::ofstream openOfstream(std::string filename){
    std::ofstream outdata;
    outdata.open(filename);
    return outdata;
  }

  void closeOfstream(std::ofstream &outdata) {
    outdata.close();
  }

  std::ostream& getCout() {
    std::ostream &out(std::cout);
    return out;
  }
}

#endif // PYROL_STREAM
