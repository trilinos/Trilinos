#ifndef PIKE_DATA_TRANSFER_HPP
#define PIKE_DATA_TRANSFER_HPP

#include "Pike_BlackBox_config.hpp"
#include <string>

namespace pike {

  class Solver;

  class DataTransfer {

  public:

    virtual ~DataTransfer() {};

    virtual std::string name() const = 0;

    virtual bool doTransfer(const pike::Solver& solver) = 0;

    virtual bool transferSucceeded() const = 0;

  };

}

#endif
