#ifndef _IFP_MATRIX_H_
#define _IFP_MATRIX_H_

class ifp_Matrix
{
public:
    virtual ~ifp_Matrix() {}

    virtual int dimrow() const = 0;
    virtual int dimcol() const = 0;

    virtual void mult(int, int, const double *, int, double *, int) const {}
    virtual void trans_mult (int, int, const double *, int, double *, int) const
      {}
};

#endif // _IFP_MATRIX_H_
