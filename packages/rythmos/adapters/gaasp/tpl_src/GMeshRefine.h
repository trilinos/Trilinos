#ifndef INC_GMeshRefine_h
#define INC_GMeshRefine_h

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
//	File:	  GMeshRefine.h
//	Class:	  GMeshRefine
//
//	Description:
//	Handles mesh refinement.
// ----------------------------------------------------------------------------
//	Author:	Jeff Sandelin, sandelin@math.colostate.edu, November, 2007
//	History:
//	<date, eg., 29May01>	<your name>, <your e-mail address>
//	<description of modifications>
// ----------------------------------------------------------------------------
#include "MemoryMgmt.h"
#include "GModelBase.h"
#include "InitGaaspOO.h"
#include <vector>
#include <ctime>

namespace GAASP {

using std::time;

class GMeshRefine
{
  public:
	//---- constructors and destructor

	GMeshRefine();
	GMeshRefine(double *, double **, TRefinementFlag, TMethodFlag, double, int, int);
	void Setup(double *, double **, TRefinementFlag, TMethodFlag, double, int, int);
	~GMeshRefine ();

	//---- functions
	void Clear ();					// "Clear" data members
	char const * const GetVersion () const		// Return version
	  { return version; }
	void Initialize();

	void refine();
	void computeATOL();
	void flagIntervals();
	void refineMesh();	
	void updateData(double **, double *, int);
	double *getNewMesh(){return newNodes;}
	int getNumberofNodes(){return nNewNodes;}

  protected:
	//---- constants

	//---- data
	bool modified;					// true if data changed
  bool isSetup_;
		
	//---- functions

  private:
	//---- constants
	static char const * const version;		// class version number

	//---- data
	int *flags;						//	Refinement flags
	bool flagsSet;					//	Are the flags set?
	int nOldNodes;					//	Number of old nodes
	int nOldIntervals;				//	Old nodes - 1
	double *oldNodes;				//	Old nodes array
	int nNewNodes;					//	Number of new nodes
	int nNewIntervals;				//	New nodes - 1;
	double *newNodes;				//	New node array
	double **intervalError;			//	Array of interval error contributions
	TMethodFlag MethodFlag;			//	dG(0) or dG(1)
	TRefinementFlag rOption;		//	Refinement option
	int sysDim;						//	System dimension
	double userTOL;					//	User defined tolerance
	double ATOL;					//	

	//---- option flags
	
	//---- vectors for initial conditions and parameters

	//---- functions
	void Copy (GMeshRefine const & object);	// Copy to this
	int binarySearch(int, int, double const);

};

} // namespace GAASP

#endif // INC_GMeshRefine_h
