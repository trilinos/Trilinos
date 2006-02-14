#ifndef SOLVER_CRD_HPP
#define SOLVER_CRD_HPP

using namespace std;

class solver_crd 
{
public: // functions
  solver_crd() { };
  virtual ~solver_crd() { };
  virtual int factor(double* A, int* rowbeg, int* colidx, int n, 
		     int scale_option=0)=0;
  virtual int solve(const int NRHS, double* RHS, double* SOL, double* TEMP)=0;
private:
protected:
  
};
#endif // SOLVER_CRD_HPP
