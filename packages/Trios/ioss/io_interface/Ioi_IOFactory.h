/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef Ioi_IOFactory_h
#define Ioi_IOFactory_h

#include <string>

namespace Ioss {
namespace Interface {
  class IOBroker;

  class IOFactory {
  public:
    virtual ~IOFactory() {};
    static IOBroker* create(const std::string& region_name);

  protected:
    explicit IOFactory();
    virtual IOBroker* make_IO(const std::string& region_name) const = 0;

  private:
    static IOFactory* factory_;
  };
}//namespace Interface
}//namespace Ioss
#endif // Ioi_IOFactory_h
