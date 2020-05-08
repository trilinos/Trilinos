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

  int _num_threads;
  int _numa;
  int _device;
  int _ngpu;

  bool _called_initialize;

public:

  // call with command line arguments
  KokkosParser(int argc, char* args[], bool print_status = false);

  // call with integer arguments
  KokkosParser(int num_threads = -1, int numa = -1, int device = -1, int ngpu = -1, bool print_status = false);

  // destructor
  ~KokkosParser() {
      // clean-up Kokkos
      if (_called_initialize) {
          this->finalize();
      } 
  };

  // initialize Kokkos if not already initialized using
  // arguments provided at object construction
  int initialize();
  
  // finalize Kokkos if this object initialized it
  // or if hard_finalize is true
  int finalize(bool hard_finalize = false);

  // retrieve arguments from a previous initialized Kokkos state
  void retrievePreviouslyInstantiatedKokkosInitArguments();

  // use this object's state to create a Kokkos object used
  // later for initialization
  Kokkos::InitArguments createInitArguments() const;

  int getNumberOfThreads() const { return _num_threads; }
  int getNuma() const { return _numa; }
  int getDeviceID() const { return _device; }
  int getNumberOfGPUs() const { return _ngpu; }

  // prints status
  void status() const;

  // prohibit using the assignment constructor
  KokkosParser& operator=( const KokkosParser& ) = delete;
  //KokkosParser( const KokkosParser& ) = delete;

};

} // Compadre

#endif
