#ifndef _IFP_GLOBALPRECON_H_
#define _IFP_GLOBALPRECON_H_

#include "ifp_Precon.h"
#include "ifp_LocalPrecon.h"

class ifp_GlobalPrecon : public ifp_Precon
{
protected:
    ifp_LocalPrecon local_precon;

public:
    ifp_GlobalPrecon() {local_precon.name = (LocalPreconName) 0;}
    virtual ~ifp_GlobalPrecon() {}

    // set the method for solving or inverting the blocks

    void localprecon(LocalPreconName b) 
        {local_precon.name = b;}
    void localprecon(LocalPreconName b, int i1) 
        {local_precon.name = b; local_precon.iarg1 = i1;}
    void localprecon(LocalPreconName b, double d1) 
        {local_precon.name = b; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, int i1, double d1) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, double d1, int i1) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, int i1, int i2) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.iarg2 = i2;}
    void localprecon(LocalPreconName b, double d1, double d2) 
      {local_precon.name = b; local_precon.darg1 = d1; local_precon.darg2 = d2;}
};

#endif // _IFP_GLOBALPRECON_H_
