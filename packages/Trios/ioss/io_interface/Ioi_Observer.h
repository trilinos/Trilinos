/* Copyright 2001 Sandia Corporation, Albuquerque, NM. */
#ifndef Ioi_Observer_h
#define Ioi_Observer_h

#include <Ioi_Events.h>

namespace sierra {
namespace Frio {
  class IOBase;
}
}

namespace Ioi {

  struct IOEvent {
    sierra::Frio::IOBase *		    m_ioBase;
    EventInterest   m_interest;
    EventState	    m_state;
  };
}

#endif // Ioi_Observer_h
