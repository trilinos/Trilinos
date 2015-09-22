#ifndef INC_GForwardSolve_h
#define INC_GForwardSolve_h

// ----------------------------------------------------------------------------
//	Developed with funding from the Department of Energy
//	Grants DE-FG02-04ER25620, DE-FG02-05ER25699, DE-FC02-07ER54909
//	and
//	National Science Foundation
//	Grants DMS-0715135, DGE-0221595003, MSPA-CSE-0434354
// ----------------------------------------------------------------------------
//	Copyright (c) 2007 Colorado State University. All rights reserved.
// ----------------------------------------------------------------------------
//	Organization:	Department of Mathematics - Estep Research Group
//					Colorado State University, Fort Collins, CO 80523 USA
//					www.math.colostate.edu/~estep/gaasp
//	Project:  Globally Accurate Adaptive Sensitivity Package
//	File:	  GForwardSolve.h
//	Class:	  GForwardSolve
//
//	Description:
//	Handles solving the forward problem for GAASP.
// ----------------------------------------------------------------------------
//	Author:	Jeff Sandelin, sandelin@math.colostate.edu, November, 2007
//	History:
//	<date, eg., 29May01>	<your name>, <your e-mail address>
//	<description of modifications>
// ----------------------------------------------------------------------------
#include "MemoryMgmt.h"
#include "GModelBase.h"
#include "InitGaaspOO.h"
#include "boost/shared_ptr.hpp"
#include <vector>

namespace GAASP {

using boost::shared_ptr;

struct TSimControl
{
	double startTime;
	double endTime;
	double stepSize;
	int nSteps;
	TMethodFlag methodFlag;
};

class GForwardSolve
{
  public:
	//---- constructors and destructor
	GForwardSolve();

	//---- Model class, TSimControl, initial conditions 
	GForwardSolve(shared_ptr<GModelBase>, TSimControl, std::vector<double>);
	//---- Model class, TSimControl, mesh, steps, initial conditions, parameters
	GForwardSolve(shared_ptr<GModelBase>, TSimControl, double *, int, std::vector<double>, std::vector<double>);
	//---- Model class, TSimControl, initial conditions, parameters, mesh 
	GForwardSolve(shared_ptr<GModelBase>, TSimControl, std::vector<double>, std::vector<double>);

  //---- Set up object 
	void Setup(shared_ptr<GModelBase>, TSimControl, std::vector<double>);
	void Setup(shared_ptr<GModelBase>, TSimControl, double *, int, std::vector<double>, std::vector<double>);
	void Setup(shared_ptr<GModelBase>, TSimControl, std::vector<double>, std::vector<double>);

	//---- Destructor
	~GForwardSolve ();

	//---- functions
	void Clear ();					// "Clear" data members
	char const * const GetVersion () const		// Return version
	  { return version; }

	//	The forward solver
	int solve();							

	double **getSolution();					//	Retrieve the solution
	double **getDerivatives();				//	Retrieve the derivatives
	int getNSteps(){return nsteps;}
	double *getMesh();						//	Retrieve the mesh
	void printSolution();					//	Write solution to file
	void printSolution(int);					//	Write solution to file
	void updateMesh(double *, int);			//	Enter a new mesh
	void updateIC(std::vector<double>);

  protected:
	//---- constants

	//---- data
	bool modified;					// true if data changed

	//---- functions
	shared_ptr<GModelBase> modelPtr;

  private:
	//---- constants
	static char const * const version;		// class version number
	int sysDim;								// system dimension
	double machEpsilon;						// machine epsilon
	int icCnt;								// IC count
	double *IC;								// Initial conditions

	//---- data
  bool isSetup_;

	//---- option flags
	TMethodFlag MethodFlag;
	
	//---- vectors for initial conditions and parameters
	std::vector<double> initialState;

	//	These can be changed to GVectors when the time comes.
	//	Note that the derivatives are only computed when the method is dG(1).
	//  For dG(0), the derivatives are zero.
	double **solution;						//	Array for solution data
	double **derivatives;					//	Array of solution derivatives
	double *tNodes;							//	Time mesh
	int nsteps;								//	Number of time steps

	void allocateSolutionSpace();

	//---- functions
	void Initialize ();						// Initialize members
	void Copy (GForwardSolve const & object);	// Copy to this

	//	Provide the solver with a guess
	void predict( double *, double *, double const, double const, int const, int const);							
	//	Newton solve
	bool newton(double *, double *, double *, double *, double const, double const, int const, int const);
	//	Forward step based on method 
	void fMethod(double *, double *, double *, double const, double const, int const);
	//	Finite difference Jacobian
	double** forwardJacobian(int const, double *, double *, double *,	double const, double const, int const);						
};

} // namespace GAASP

#endif // INC_GForwardSolve_h
