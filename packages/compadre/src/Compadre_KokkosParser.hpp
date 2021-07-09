#ifndef _COMPADRE_KOKKOSPARSER_HPP_
#define _COMPADRE_KOKKOSPARSER_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"

namespace Compadre {

/*! \class KokkosParser
    \brief Class handling Kokkos command line arguments and returning parameters.
*/
class KokkosParser {

private:

  bool _called_initialize;

  // prevent default constructor
  KokkosParser();

public:

  // call with command line arguments
  KokkosParser(int argc, char* args[], bool print_status = false);

  // call with std::vector of std::string's
  KokkosParser(std::vector<std::string> args, bool print_status = false);

  // call for default arguments
  KokkosParser(bool print_status = false);

  // destructor
  ~KokkosParser() {
      // clean-up Kokkos
      if (_called_initialize) {
          this->finalize();
      } 
  };

  // initialize Kokkos if not already initialized using
  // arguments provided at object construction
  int initialize(int argc, char*[], bool print_status = false);
  
  // finalize Kokkos if this object initialized it
  // or if hard_finalize is true
  int finalize(bool hard_finalize = false);

  // prints Kokkos configuration
  void status() const;

  // prohibit using the assignment constructor
  KokkosParser& operator=( const KokkosParser& ) = delete;

};

} // Compadre

#endif
