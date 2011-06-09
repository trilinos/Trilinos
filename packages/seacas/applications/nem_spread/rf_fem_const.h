/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*
*
*    Definition of globally occurring constants relating
*    to the details of the FEM problem.
*/

/* Geometric parameters */

#define CARTESIAN	1
#define CYLINDRICAL	2
#define SPHERICAL	3

/*
 * Variable Types : Each unknown in salsa is identified by a variable
 *                  type and variable subtype pair.  This is the
 *                  official list.
 *
 *     Future plans: Make this into an enumerated type
 */

#define VELOCITY1	   0
#define VELOCITY2	   1
#define VELOCITY3	   2
#define PRESSURE	   3
#define TEMPERATURE	   4
#define MASS_FRACTION	   5     /* subtype is the gas species number    */
#define SURFACE   	   6
#define MESH_DISPLACEMENT1 7
#define MESH_DISPLACEMENT2 8
#define MESH_DISPLACEMENT3 9
#define SITE_FRACTION      10    /* subtype is the surface species num   */
#define SITE_DENSITY       11    /* subtype is the surf phase num        */
#define BULK_MOLE_FRACTION 12    /* subtype is the bulk species num      */
#define FLUX               10
/*
 * NORMAL and TANGENT velocity variables must be distinct from
 * standard velocity variables so that BC exact functions can
 * distinguish them from variables for VELOCITY1, etc.
 */
#define NORMAL_VELOCITY    20
#define TANGENT_VELOCITY1  21
#define TANGENT_VELOCITY2  22
/*
 * NORMAL and TANGENT VELOCITY_EQNs must have same definitions
 * as VELOCITY1, VELOCITY2, and VELOCITY3 (although not necessarily
 * in the same order).  These defined names are used to determine
 * which momentum equation to fill with normal and tangent velocity
 * boundary conditions.
 */
#define NORMAL_VELOCITY_EQN   0
#define TANGENT_VELOCITY1_EQN 1
#define TANGENT_VELOCITY2_EQN 2

/*
 * Constant based on the number of entries above:
 *
 *   (note: this is stuck at 7 until we upgrade the mapping arrays to
 *          handle types 7 to 13)
 */

#define MAX_VARIABLE_TYPES 7


/* Type of velocity-pressure elements */

#define	 Q2Q1  10 /* Quadratic velocity, linear pressure; LBB stable     */
#define	 Q2Q2  15 /* Quadratic velocity, Quadratic pressure; LBB unstable*/
#define  pQ2Q1 20 /* Pseudo version of Q2Q1; superelement; LBB stable    */
#define  pQ2P0 30 /* Piecewise constant pressure over superelement;      */
#define  Q1Q1  40 /* Bilinear velocity and pressure; LBB unstable        */
#define  Q1P0  50 /* Bilinear vel., piecewise const. press.;LBB unstable */

/*
 *    Selection of Problem and Mechanism
 *        - Valid values for the global variable, ProblemType
 */

#define  ENERGY_DIFF            0    /* Conduction heat transfer          */
#define  ENERGY_CONV_DIFF       1    /* Conduction with convection        */
#define  MASS_DIFF              2    /* Mass diffusion transport          */
#define  MASS_CONV_DIFF         3    /* Mass diffusion with convection    */
#define  ENERGY_MASS_DIFF       4    /* Conduction and Mass diff.         */
#define  ENERGY_MASS_CONV_DIFF  5    /* Conduction and Mass diff. w/conv. */
#define  STOKES_FLOW            6    /* No inertial terms in momentum     */
#define  FLUID_FLOW             7    /* Solve fluids problem only         */
#define  FLUID_FLOW_ENERGY      8    /* Solve fluids problem and heat tran*/
#define  FLUID_FLOW_MASS        9    /* Solve fluids problem and Mass tran*/
#define  TURB_FLOW              10   /* Solve turbulent flow only         */
#define  TURB_FLOW_ENERGY       11   /* Solve turbulent flow and heat tran*/
#define  TURB_FLOW_MASS         12   /* Solve turbulent flow and Mass tran*/
#define  WHOLE_BANANA           13   /* Every thing!                      */

/*
 *    Selection of linearity choices
 *       - Valid values for the global variable, ProblemLinearity
 */

#define  DEFAULT       0    /* Solve problem with default linearity defined   */
                            /* by the ProblemType                             */
#define  LINEAR        1    /* Force problem to be solved as a linear problem */
#define  NONLINEAR     2    /* Force problem to be solved as a nonlinear prb. */

/*
 *  Selection of Time Integration Method
 *       Valid values for the global variable, TimeIntegration
 */

#define  STEADY          0    /* Steady state solution method               */
#define  TRANSIENT       1    /* Accurate transient solution method         */
#define  PSEUDO          2    /* Pseudo transient method (not time accurate)*/
#define  CONTINUATION    3    /* Continuation of steady state solutions     */
#define  OPTIMIZATION    4    /* Optimization of steady state solutions     */

/*
 * Selection of the strategy for selection of the time step control
 *       Valid values for the global variable, TimeStepControl
 */

#define  VARIABLE        100  /* Dynamic time step control strategy         */
#define  CONSTANT        101  /* Constant time increments                   */


/* Selection of Multicomponent transport algorithm */

#define  MIXTURE_AVG     1    /* Mixture averaged approximation             */
#define  STEFAN_MAXWELL  2    /* Oran-Boris approx. to Stefan-Maxwell form. */
#define  DIXON_LEWIS     3    /* Full multicomponent formulation            */

/*
 *  Equation of State
 */

#define BoussinesQ      101


/* Temp defs */

#define LINEAR          1
#define QUADRATIC       2

/******************************************************************************/
/*		EXTERN DECLARATIONS FOR GLOBALS DEFINED IN rf_fem.h	      */
/*					VERSION 1.22			      */
/******************************************************************************/

/* Geometric Paramters  */
extern
int  CoordinateSystem;  /* Indicates type of coordinate system
			   (see fem_const.h)*/

/* FEM Interpolation Parameters  */

extern
int  VelocityPressure;  /* Indicates which element type is used             */
extern	                /* for velocity and pressure interpolation.         */
int  Velocity;		/* Indicates which type of interpolation is used    */
extern			/* for velocity. (set by value of VelocityPressure) */
int  Pressure;		/* Indicates which type of interpolation is used    */
extern			/* for pressure. (set by value of VelocityPressure) */
int  Temperature;       /* Indicates which type of interpolation is used    */
extern			/* for temperature.                                 */
int  MassFraction;      /* Indicates which type of interpolation is used    */
                        /* for mass fraction and density.                   */

/* Flags to indicate the problem type */

extern
int  FluidMechanics,	/* If TRUE, solve momentum transport and continuity  */
     MassTransfer,	/* If TRUE, solve species transport equation         */
     HeatTransfer,	/* If TRUE, solve thermal energy transport equation  */
     Turbulence;        /* If TRUE, include turbulence transport in some form*/

extern
int  VelocitySet;	/* If TRUE, there is a prescribed velocity field
				    and therefore no momentum equation needs
				    to be solved 	                    */
extern
int  TemperatureSet;	/* If TRUE, there is a rescribed temperature field
				    and therefore, no energy equations
				    needs to be solved                      */
extern
int  Incompressible,	/* If TRUE, no mixture density variation            */
     InviscidFluid,	/* If TRUE, there is no diffusive transport of
				    energy,  momentum, or species   	    */
     MomentumAdvec,	/* If TRUE, there is advection of momentum
			            in the momentum transport eq.   	    */
     EnergyAdvec,	/* IF TRUE, there is advection of energy in the
			            energy equation	                    */
     SpeciesAdvec;	/* If TRUE, there is advection of species in the
			            species conservation equation           */
extern
int  EnergyRx,		/* If TRUE include energy reaction source terms     */
     SpeciesRx,		/* If TRUE include species reaction source terms    */
     ThermalDiffusion,	/* If TRUE include thermal diffusion effects        */
     ChemicalEq,	/* If TRUE thermodynamic chemical equilibrium exists*/
     DiluteApprox;	/* If TRUE use dilute approx. to mass transport     */

extern
int  NonLinear;		/* if TRUE use nonlinear solution methods           */

/* Parameters to select Problem type  */

extern
int  ProblemType,	/* Select type of problem to be solved		    */
     ProblemLinearity,  /* Select either default, linear or nonlinear       */
                        /* coupling of dependent variables                  */
     Multicomponent;	/* Select fomulation for multicomponent transport   */


/* Flags to select matrix terms for discrete transport equations */

extern
int MomentumMassMatrix,	/* if TRUE the individual terms for these         */
    PressureMassMatrix,	/* matrices contribute to the global coefficient  */
    SpeciesMassMatrix,	/* matrix.					  */
    EnergyMassMatrix,
    MomemtumAdvecMatrix,
    EnergyAdvecMatrix,
    SpeciesAdvecMatrix,
    MomentumDiffMatrix,
    EnergyDiffMatrix,
    SpeciesDiffMatrix,
    SpeciesThermalDiffMatrix,
    DivergenceMatrix,
    MomentumBodyForceMatrix;


/*
 *
 * Flags to enable different stabilization formulations.
 *
 * NOTE: Petrov-Galerkin Pressure stabilization is always on.
 *
 */

extern
int SUPG_Flag; /* If true Strealine Upwind Petrov-Galerkin convection
                  stabilization is enabled. Default = FALSE */


/* Boundary Condition information */

extern
int Num_BC; 	      /* Number of boundary conditions defined            */
extern
int Num_GS;           /* Number of generalized surface definitions        */
extern
int Num_Fn_Data;      /* Number of user functions in the input file       */
extern
int Num_Enclosures;   /* Number of radiation enclosures                   */
extern
int Num_Tfunc;
extern
int Num_Vfunc;

/* Paramters to select time integration technique                           */

extern
int  TimeIntegration,     /* Select time integration method                   */
     TimeIntOrder,        /* Select order of time integration method          */
     TimeStepControl;     /* Flag to enable dynamic time step control         */
extern
int  MaxTimeSteps;	  /* Maximum number of time steps in transient soution*/
extern
double Delta_t0,	  /* Initial time step for transient solution         */
       TimeMax,		  /* Maximum upper limit on time in transient solution*/
       ParamInit;  	  /* Inital value of the parameter or a flag value    */

/* Internal Flags that indicate what needs to be calculated within the
   Residual and Coefficient Matrix Function                                 */

   /* The following flags can have TRUE and FALSE values */
extern
int  NumRegUnk_PN ;      /* Number of regular unknowns defined at each      */
                         /* global node.  This is defined as the number of  */
                         /* unknowns which have basis functions defined at  */
                         /* all global nodes in the element.  Unknowns      */
                         /* whose basis function interpolations are really  */
                         /* subparametrizations of the element are not      */
                         /* included, here.                                 */
extern
int  SurfInt,           /* Indicates that a surface integral calculation is */
                        /* necessary for some equation for some elements    */
     SurfIntMomentum,   /* Indicates that a surface integral calculation is */
                        /* necessary for the momentum equation (all of them)*/
                        /* for some elements.                               */
     SurfIntEnergy,     /* Indicates that a surface integral calculation is */
                        /* necessary for the energy equation                */
                        /* for some elements.                               */
     SurfIntSpecies;    /* Indicates that a surface integral calculation is */
                        /* necessary for all of the species equations       */
                        /* for some elements.                               */


/* Internal Parameters that calculate information needed in the
   Residual and Coefficient Matrix Function                                 */
extern
int  NumSolnVelocities; /* Number of velocities to be solved for            */
                        /* (0, 1, 2, or 3)                                  */
extern
int SolveThetaVelocity; /* if TRUE Solve for the theta component of velocity */
			/* in an axisymmetric cylindrical coordinate system  */


/* Information on the number of unknowns (variables) which are define       */

extern
int GNumUnknowns;       /* Total Number of unknown variables in problem  */
extern
int NumUnknowns;        /* Number of unknown variables updated by this   */
extern			/* processor                                     */
int NumExtUnknowns;     /* Number of external variables which are        */
extern                  /* information computed on this proc in residual */
int NumIntUnknowns;     /* Number of internal variables which use only   */
extern			/* copied and stored on the local processor      */
int MaxVarPerNode;      /* Maximum number of unknowns at any node on the */
extern			/* processor                                     */
unsigned char *Num_Unknowns_Node;
			/* The number of unknowns which are define at    */
extern			/* node                                          */
int Num_Var_In_Type[MAX_VARIABLE_TYPES];
			/* The number of variables of this type. (i.e.   */
			/* for the MASS_FRACTION type there are KkGas    */
extern			/* species.					 */
int *First_Unknown;	/* Index to start of the first unknown at each   */
			/* node. To determine where an unknown ar node i */
			/* exists in the solution vector, use:           */
			/* First_Unknown[i] and Variable_Mask[i].        */
extern
signed char *First_Y;        /* Incremental index to start of first mass      *
			      * fraction unknown at node i.                   *
			      * (returns -1 if there is no Y at i )           */

extern
int *Index_P;	        /* Index location in solution vector for the     */
			/* pressure variable which is defined at this    */
			/* local node. If pressure is not defined ->     */
			/* value= -1                                     */
extern
double Ptherm;          /* The thermodynamic pressure for the problem.   */

extern
double RU_Const;        /* The universal gas constant (cal/mole*K).      */

/******************************************************************************/
/*		PROTOTYPE DEFINITIONS FOR GLOBAL FEM FUNCTIONS 		      */
/******************************************************************************/

/*--------------- rf_fanin.c -----------------*/

extern int gather_nodal_var (
                        int    num_nodes,
       	  		double sol_vec[],
         		int    gnodes[],
         		int    var_no,
         		int    k,
         		int   *N_proc,
         		int   *gindex[],
         		float *gvec[]
		  );

extern int gather_elem_var (
			    double **sol_vec,
			    int    gelems[],
			    int    k,
			    int    *N_proc,
			    int    gindex[],
			    float  gvec[]
			    );


/******************************************************************************/
