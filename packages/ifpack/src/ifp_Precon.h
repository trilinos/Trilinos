#ifndef _IFP_PRECON_H_
#define _IFP_PRECON_H_

class ifp_Precon
{
public:
    virtual ~ifp_Precon() {}

    virtual void apply  (int, int, const double *, int, double *, int) {}
    virtual void applyt (int, int, const double *, int, double *, int) {}
    virtual void applyr (int, int, const double *, int, double *, int) {}
    virtual void applyrt(int, int, const double *, int, double *, int) {}
    virtual void applyl (int, int, const double *, int, double *, int) {}
    virtual void applylt(int, int, const double *, int, double *, int) {}
    virtual void applyc (int, int, const double *, int, double *, int) {}
    virtual void applyct(int, int, const double *, int, double *, int) {}
};

#endif // _IFP_PRECON_H_
