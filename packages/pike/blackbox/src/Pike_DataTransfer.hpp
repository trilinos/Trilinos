#ifndef PIKE_DATA_TRANSFER_HPP
#define PIKE_DATA_TRANSFER_HPP

#include "Pike_BlackBox_config.hpp"
#include <string>
#include <vector>

namespace pike {

  class Solver;
  class BlackBoxModelEvaluator;

  class DataTransfer {

  public:

    virtual ~DataTransfer() {};

    virtual std::string name() const = 0;

    virtual bool doTransfer(const pike::Solver& solver) = 0;

    virtual bool transferSucceeded() const = 0;

    virtual const std::vector<std::string>& getSourceModelNames() const = 0;
    
    virtual const std::vector<std::string>& getTargetModelNames() const = 0;
  };

}

#endif
