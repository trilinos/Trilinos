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
/******************************************************************************/
/*            Prototypes for user defined functions                           */
/******************************************************************************/

extern void eval_vfunc (
     double     *value,          /* Value to be returned -Why not just
                                    declare this a double function?          */
     double       temp,          /* value to be interpolated in the
                                    lookup table                             */
     int           id            /* integer id of the lookup table           */
);

extern void eval_tfunc(double *value, double time, int id);

extern void user_out (
     int         status,          /* Flag indicating spot in program where
                                     subroutine is being called               */
     double     *soln_vec,        /* Solution vector                          */
     double      time,            /* Time                                     */
     int         time_step_num    /* Time step number, if applicable          */
);


extern void user_density(

     double      *rho,             /* density at (x,y,z)        */
     double       temperature,     /* Temperature at (x,y,z)    */
     double       X_k[],           /* Vector of  mole fraction
                                      variables at the current position
                                      X_k[Num_species]          */
     double       Ptherm,          /*  Thermodynamic pressure   */
     double       x,
     double       y,
     double       z,               /* location (x,y,z) */
     MATSTRUCT_PTR matID_ptr       /* material ID number        */
);


extern void user_visc (
     double      *viscosity,       /* viscosity at (x,y,z)      */
     double       temperature,     /* Temperature at (x,y,z)    */
     double       X_k[],           /* Vector of  mole fraction
                                       variables at current position
                                       X_k[Num_species]         */
     double       Ptherm,          /*  Thermodynamic pressure   */
     double       x,
     double       y,
     double       z,               /* location (x,y,z) */
     MATSTRUCT_PTR matID_ptr      /* Pointer to the current material
                                      structure                   */
);

extern void user_Cp (
     double      *Cp,              /* specific heat at (x,y,z)  */
     double       temperature,     /* Temperature at (x,y,z)    */
     double       X_k[],           /* Vector of  mole fraction
                                       variables at current position
                                       X_k[Num_species]         */
     double       Ptherm,          /* Thermodynamic pressure   */
     double       x,
     double       y,
     double       z,               /* location (x,y,z) */
     MATSTRUCT_PTR matID_ptr      /* Pointer to the current material
                                      structure                   */
);

extern void user_cond (
     double      *conductivity,    /* thermal conductivity
                                               at (x,y,z)                */
     double       temperature,     /* Temperature at (x,y,z)    */
     double       X_k[],           /* Vector of mole fractions
                                      variables at current position.
                                      X_k[Num_species] */
     double       Ptherm,          /* Thermodynamic pressure   */
     double       x,
     double       y,
     double       z,               /* location (x,y,z) */
     MATSTRUCT_PTR matID_ptr      /* Pointer to the current material
                                      structure                   */
);

extern void user_hcoeff (double *hcoeff, double temperature, double tref,
                  int ss_id, double time);

extern void user_tref (double *tref, double temperature, int ss_id,double time);

extern void user_formf (double *formf,double temperature,int ss_id,double time);

extern void user_continuation (
     double      *ContinuationParam  /* Address of Continuation Parameter */
);

#ifdef OPTIMIZER
extern int user_optimization (
     double      *optimizationParam  /* Address of Continuation Parameter */
);

extern double obj_func_eval (
     double opt_index,      /* Optimization step number */
     double status,         /* status -- see user_out */
     double depo_rate,      /* Average deposition rate from f_xy_spin_average */
     double uniform_area    /* Area of good crystal from f_xy_spin_average */
);
#endif

/************************* END of rf_user.h ***********************************/
