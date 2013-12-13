#ifndef PIKE_DATA_TRANSFER_HPP
#define PIKE_DATA_TRANSFER_HPP

namespace pike {

  class DataTransfer {

  public:

    virtual ~DataTransfer() {};

    virtual void doTransfer() = 0;

    virtual bool transferSucceeded() = 0;

  };

}

#endif
