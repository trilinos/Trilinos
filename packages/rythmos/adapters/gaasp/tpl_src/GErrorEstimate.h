#ifndef INC_GErrorEstimate_h
#define INC_GErrorEstimate_h

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
//	File:	  GErrorEstimate.h
//	Class:	  GErrorEstimate
//
//	Description:
//	Handles the error estimates for GAASP
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

class GErrorEstimate
{
  public:
	//---- constructors and destructor

	//---- 
	GErrorEstimate();
	GErrorEstimate(double **, double **, double **, double *, int, int, shared_ptr<GModelBase>, TMethodFlag);
	void Setup(double **, double **, double **, double *, int, int, shared_ptr<GModelBase>, TMethodFlag);
	~GErrorEstimate ();

	//---- functions
	void Clear ();					// "Clear" data members
	char const * const GetVersion () const		// Return version
	  { return version; }
	shared_ptr<GModelBase> modelPtr;

	//	Error estimates
	int compute();							
	void updateForwardSolution();
	void updateForwardDerivatives();
	void updateAdjointSolution();
	void updateMesh();
	void updateData(double *, int, double **, double **, double **);
	double getErrorEstimate();
	double **getIntervalContributions();
	void printErrorEstimates();

  protected:
	//---- constants

	//---- data
	bool modified;					// true if data changed
  bool isSetup_;

	//---- functions

  private:
	//---- constants
	static char const * const version;		// class version number
	int sysDim;								// system dimension

	//---- data
	int static const gaussRule = 5;
	double t[5];
	double w[5];
	double tau[2];	//	Radau points
	double dt;		//	Delta t
	double dtau;	//	delta tau

	//---- option flags
	TMethodFlag MethodFlag;
	
	//---- vectors for initial conditions and parameters
	std::vector<double> initialState;

	//	These can be changed to GVectors when the time comes.
	int nsteps;								//	Number of time steps
	int ssteps;								//	Number of spline steps (determined by adjoint)
	double *tNodes;							//	Time mesh
	double *sNodes;							//	Adjoint mesh
	double **fsoln;							//	Array for forward solution data
	double **sder;							//	Array for forward solution derivatives
	double **dsoln;							//	Array for adjoint solution data
	double **intervalError;					//	Error by equation and time step
	double *intervalContributions;			//	Interval  contribution
	double *discError;						//	Interval discretization error
	double *quadError;						//	Interval quadrature error

	double **dualy2;						//	Second erivatives for the spline

	void allocateSolutionSpace();

	//---- functions
	void Initialize ();						// Initialize members
	void Copy (GErrorEstimate const & object);	// Copy to this
	void generateAdjointMesh();

	//---- functions used by the error computation
	void discretizationError(int,double,double,double);
	void quadratureError(int,double,double,double);

	double *interpolant( int, int, int, int );
	double gamma( int, double );
	void yrad( int, int, double * );
	void ybar( int, double, double * );
	void fbar( int, double, double * );
	void fBAR( int, double, double * );

	friend void spline (double *, double *, int, double * const);
	friend int splint (double  const * const, double  const * const, double  const * const,
		int  const, double  const, double * const);
	friend double projection(double, double, double, double, double);

};

} // namespace GAASP

#endif // INC_GErrorEstimate_h
