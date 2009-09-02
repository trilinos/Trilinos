#ifndef INC_GAdjointSolve_h
#define INC_GAdjointSolve_h

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
//	File:	  GAdjointSolve.h
//	Class:	  GAdjointSolve
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

struct TAdjointControl
{
	const double startTime;		// Simulation start time
	const double endTime;		// Simulation end time
	const double stepSize;		// Time step size
	const int nSteps;			// Number of steps
	TErrorFlag errorFlag;		// Quantity of interest
	TMethodFlag methodFlag;		// Low/High order method
	const int errcomp;			// Quantity of interest
};


class GAdjointSolve
{
  public:
	//---- constructors and destructor
	GAdjointSolve ();	
	//---- function, method, forward solution, mesh
	GAdjointSolve (shared_ptr<GModelBase>, TAdjointControl, double **, double *);	
	//---- function, method, forward solution filename
	GAdjointSolve (shared_ptr<GModelBase>, TAdjointControl, std::string const &);

	void Setup (shared_ptr<GModelBase>, TAdjointControl, double **, double *);	
	void Setup (shared_ptr<GModelBase>, TAdjointControl, std::string const &);

	~GAdjointSolve ();

	//---- data

	double **solution;						//	Array for adjoint solution data
  bool isSetup_;

	//---- functions
	shared_ptr<GModelBase> modelPtr;

	void Clear ();					// "Clear" data members
	char const * const GetVersion () const		// Return version
	  { return version; }

	//---- The adjoint solver
	int solve();	
	void printSolution();
	void printSolution(int);
	double **getSolution();					//	Retrieve the solution
	void updateData(double *, int, double **);

protected:
	//---- constants

	//---- data
	bool modified;					// true if data changed

	//---- functions

  private:
	//---- constants
	static char const * const version;		// class version number
	int sysDim;								// system dimension
	double machEpsilon;						// machine epsilon

	//---- data
	double **spln;							// Forward solution spline

	//---- spline variables for Becky's spline.
//	double **ypp;							// second derivative
//	double ypval;							// first derivative at a specific time
//	double yppval;							// second derivative at a specific time
	//---- spline variables for Jeff's spline
	double **y2;							// second derivative

	//---- option flags
	TMethodFlag MethodFlag;
	TErrorFlag ErrorFlag;
	
	//---- vectors for initial conditions and parameters
	std::vector<double> initialState;

	//	These can be changed to GVectors when the time comes.
	double *tNodes;							//	Forward solution time mesh
	double *dNodes;							//	Adjoint solution time mesh
	double **fsoln;							//	Forward solution
	//	Note:	Not sure why we need the forward solution when we have the spline data.
	//			Checking into this.

	int nsteps;								//	Number of forward solution time steps
	int ssteps;								//	Total steps including intermediate
	int errComp;							//	Specific component of a system
	double endTime;							//  End time of the forward approximation

	int putMesh();							//	Add a mesh to the class
	void allocateSolutionSpace();

	//---- functions
	void Initialize ();						// Initialize members
	void Copy (GAdjointSolve const & object);	// Copy to this
	double *psi(int const);
	void generateAdjointMesh();

	void calcDerivs();						//	Calculate the forward approximation derivatives
	double *psi();

	double **fdMethod( double *, double *, double tnew, double, double, int, int);
	double **adjointJacobian(int const, double *, double *, double *,	double const, double const, int const);

	friend void spline (double *, double *, int, double * const);
	friend int splint (double  const * const, double  const * const, double  const * const,
		int  const, double  const, double * const);

};

} // namespace GAASP

#endif // INC_GAdjointSolve_h
