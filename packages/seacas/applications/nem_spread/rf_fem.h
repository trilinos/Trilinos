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
*	Include file for globally occuring parameters and flags
*	specific to the FEM problem.
*/

/* Geometric Paramters  */

int  CoordinateSystem;  /* Indicates type of coordinate system
			   (see fem_const.h)*/

/* FEM Interpolation Parameters (see fem_const.h) */

int  VelocityPressure;  /* Indicates which element type is used             */
	                /* for velocity and pressure interpolation.         */
int  Velocity;		/* Indicates which type of interpolation is used    */
			/* for velocity. (set by value of VelocityPressure) */
int  Pressure;		/* Indicates which type of interpolation is used    */
			/* for pressure. (set by value of VelocityPressure) */
int  Temperature;       /* Indicates which type of interpolation is used    */
			/* for temperature.                                 */
int  MassFraction;      /* Indicates which type of interpolation is used    */
                        /* for mass fraction and density.                   */

/* Flags to select problem type */

int  FluidMechanics,	/* If TRUE, solve momentum transport and continuity  */
     MassTransfer,	/* If TRUE, solve species transport equation         */
     HeatTransfer,      /* If TRUE, solve thermal energy transport equation  */
     Turbulence;        /* If TRUE, include turbulence transport in some form*/


int  VelocitySet;	/* If TRUE, there is a rescribed velocity field
				    and therefore no momentum equation needs
				    to be solved 	                    */
int  TemperatureSet;	/* If TRUE, there is a rescribed temperature field
				    and therefore, no energy equations
				    needs to be solved                      */

int  Incompressible,	/* If TRUE, no mixture density variation            */
     InviscidFluid,	/* If TRUE, there is no diffusive transport of
				    energy,  momentum, or species   	    */
     MomentumAdvec,	/* If TRUE, there is advection of momentum
			            in the momentum transport eq.   	    */
     EnergyAdvec,	/* IF TRUE, there is advection of energy in the
			            energy equation	                    */
     SpeciesAdvec;	/* If TRUE, there is advection of species in the
			            species conservation equation           */

int  EnergyRx,		/* If TRUE include energy reaction source terms     */
     SpeciesRx,		/* If TRUE include species reaction source terms    */
     ThermalDiffusion,	/* If TRUE include thermal diffusion effects        */
     ChemicalEq,	/* If TRUE thermodynamic chemical equilibrium exists*/
     DiluteApprox;	/* If TRUE use dilute approx. to mass transport     */

int  NonLinear;		/* if TRUE use nonlinear solution methods           */

/* Parameters to select Problem type (see fem_const.h)*/

int  ProblemType,	/* Select type of problem to be solved		    */
     ProblemLinearity,  /* Select either default, linear or nonlinear       */
                        /* coupling of dependent variables                  */
     Multicomponent;	/* Select fomulation for multicomponent transport   */


/* Flags to select matrix terms for discrete transport equations */

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

int SUPG_Flag = FALSE; /* If true Strealine Upwind Petrov-Galerkin
                          convection stabilization is enabled */

/*
 *           Boundary Condition information
 */

int Num_BC;		 /* Number of boundary conditions defined            */
int Num_GS;              /* Number of generalized surface definitions        */
int Num_Enclosures;      /* Number of radiation enclosures                   */
int Num_Fn_Data;         /* Number of user functions in the input file       */
int Num_Tfunc;
int Num_Vfunc;


/*
 *           Paramters to select time integration technique
 */

int  TimeIntegration,     /* Select time integration method                   */
     TimeIntOrder,        /* Select order of time integration method          */
     TimeStepControl;     /* Flag to enable dynamic time step control         */

int  MaxTimeSteps;	  /* Maximum number of time steps in transient soution*/

double Delta_t0,	  /* Initial time step for transient solution         */
       TimeMax,		  /* Maximum upper limit on time in transient solution*/
       ParamInit;	  /* Initial value of the parameter of a flag value   */

/*
 *           Internal Flags that indicate what needs to be calculated
 *           within the Residual and Coefficient Matrix Function
 */

             /* The following flags can have TRUE and FALSE values */

int  NumRegUnk_PN ;        /* Number of regular unknowns defined at each      */
                           /* global node.  This is defined as the number of  */
                           /* unknowns which have basis functions defined at  */
                           /* all global nodes in the element.  Unknowns      */
                           /* whose basis function interpolations are really  */
                           /* subparametrizations of the element are not      */
                           /* included, here.                                 */

int  SurfInt,             /* Indicates that a surface integral calculation is */
                          /* necessary for some equation for some elements    */
     SurfIntMomentum,     /* Indicates that a surface integral calculation is */
                          /* necessary for the momentum equation (all of them)*/
                          /* for some elements.                               */
     SurfIntEnergy,       /* Indicates that a surface integral calculation is */
                          /* necessary for the energy equation                */
                          /* for some elements.                               */
     SurfIntSpecies;      /* Indicates that a surface integral calculation is */
                          /* necessary for all of the species equations       */
                          /* for some elements.                               */

/*
 *            Internal Parameters that calculate information needed in the
 *            Residual and Coefficient Matrix Function
 */

int  NumSolnVelocities;  /* Number of velocities to be solved for             */
                         /* (0, 1, 2, or 3)                                   */

int SolveThetaVelocity;  /* if TRUE Solve for the theta component of velocity */
			 /* in an axisymmetric cylindrical coordinate system  */


/*
 *             Information on the number of unknowns (variables)
 *             which are defined in the problem
 */

int GNumUnknowns;            /* Total Number of unknown variables in problem  */
int NumUnknowns;             /* Number of unknown variables updated by this   */
			     /* processor                                     */
int NumIntUnknowns;          /* Number of internal variables which only use   */
			     /* unknowns owned by this proc to calculate      */			     /* the residual vector			      */
int NumExtUnknowns;          /* Number of external variables which are        */
			     /* copied and stored on the local processor      */
int MaxVarPerNode;           /* Maximum number of unknowns at any node on the */
			     /* processor                                     */
unsigned char *Num_Unknowns_Node = NULL;
			     /* The number of unknowns which are define at    */
			     /* node                                          */
int Num_Var_In_Type[MAX_VARIABLE_TYPES];
			     /* The number of variables of this type. (i.e.   */
			     /* for the MASS_FRACTION type there are KkGas    */
			     /* species.				      */
int *First_Unknown = NULL;   /* Index to start of the first unknown at each   */
			     /* node. To determine where an unknown ar node i */
		 	     /* exists in the solution vector, use:           */
			     /* First_Unknown[i] and Variable_Mask[i].        */

signed char *First_Y = NULL; /* Incremental index to start of first mass      *
			      * fraction unknown at node i.                   *
			      * (returns -1 if there is no Y at i )           */

int *Index_P = NULL;         /* Index location in solution vector for the     */
			     /* pressure variable which is defined at this    */
			     /* local node. If pressure is not defined ->     */
			     /* value= -1                                     */

double Ptherm;               /* The thermodynamic pressure for the problem.   */

double RU_Const;             /* The universal gas constant (cal/mole*K).      */
