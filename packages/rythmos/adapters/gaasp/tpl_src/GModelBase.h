#ifndef INC_GModelBase_h
#define INC_GModelBase_h

#include <vector>

namespace GAASP {

class GModelBase
{
public:
	// Functions
  virtual ~GModelBase() {};
	virtual int calcDerivs(double *yin, double *yout, double t) = 0;
	virtual void updateParameters(std::vector<double>) = 0;
	virtual int getDim() = 0;

};
 
} // namespace GAASP
#endif
