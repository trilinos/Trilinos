#ifndef INC_initGaaspOO_h
#define INC_initGaaspOO_h

namespace GAASP {

enum TSolveFlag
{
    Solve_Forward,	// 0  Solve forward only
    Solve_Adjoint,	// 1  Solve adjoint and compute error
    Solve_Error		// 2  Solve forward, adjoint, and compute error
};

enum TMethodFlag
{
    Method_DG0,		// 0 dG(0)/cG(1)
    Method_DG1		// 1 dG(1)/cG(2)
};

enum TRefinementFlag
{
    Refine_None,			// 0  No Refinement
    Refine_ByContribution,	// 1  Refine if interval contrib > tolerance
    Refine_Probabilistic1,	// 2  By contribution
    Refine_Probabilistic2,	// 3  By time and contribution
    Refine_WeightedAdaptive	// 4  Weighted adapative strategy
};

enum TErrorFlag
{
    Error_End,		// 0  Error at end time
    Error_Avg,		// 1  Average error over time
	Error_Dis		// 2  Average error weighted by distance (Lorenz)
};

// Flags for mesh types.
enum TMeshFlag
{
    Mesh_None, 
    Mesh_Forward, 
    Mesh_Adjoint, 
    Mesh_CG1, 
    Mesh_CG2
};

enum TPrintFlag		// Should get rid of this when odeSolver is class !!!Becky
{
	Print_None,				// 0
	Print_Solution,			// 1
	Print_SolutionAndError,	// 2
	Print_All				// 3  Includes interval contributions
};
/*
struct TOptionFlags
{
    const TSolveFlag solveFlag;
    const TMethodFlag methodFlag;
    const TRefinementFlag refineFlag;
    const TErrorFlag errorFlag;
    TPrintFlag printFlag;		// This may change - reason given below

    TOptionFlags (	// Set default values for option flags
	const TSolveFlag useSolveFlag,
	const TMethodFlag useMethodFlag,
	const TRefinementFlag useRefineFlag,
	const TErrorFlag useErrorFlag,
	TPrintFlag usePrintFlag )	// This is correlated with solveFlag
	: solveFlag (useSolveFlag),
	  methodFlag (useMethodFlag),
	  refineFlag (useRefineFlag),
	  errorFlag (useErrorFlag),
	  printFlag (usePrintFlag) 
    {}
};

struct TSimulationControl
{
    const double startTime;
    const double endTime;
    const double timeStepSize;
    const double tolerance;
    const unsigned short maximumRefinementCycles;
    const unsigned short systemDimension;
    const unsigned short component;	// Remove?  !!!Becky 

    TSimulationControl (	// Default values - will result in no simulation
	const double useStartTime,
	const double useEndTime,
	const double useTimeStepSize,
	const double useTolerance,
	const unsigned short useMaxRefinementCycles,
	const unsigned short useSystemDimension,
	const unsigned short useComponent )
	: startTime (useStartTime),
	  endTime (useEndTime),
	  timeStepSize (useTimeStepSize),
	  tolerance (useTolerance),
          maximumRefinementCycles (useMaxRefinementCycles),
	  systemDimension (useSystemDimension),
          component (useComponent)
    {}
};
*/

} // namespace GAASP

#endif // INC_initGaaspOO_h
