// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <math.h>
#include <stdlib.h>
#include "calc_decomp_cuts.h"

#include <stdio.h>

/* ========================================================================== */


/* Note, the TOTAL amount of data communicated (summed over all
   processes and all INTERNAL boundaries) is computed and stored in
   the variable I, the minimum of which is I0. */

namespace PAMGEN_NEVADA {


long long dom_decomp_2d(const long long Nx, const long long Ny,
		  const long long Np, long long *pNGx, long long *pNGy){

  long long rx, ry, I, rxs, rx_min, rx_max;
  long long rx0 = 1;
  long long ry0 = 1;
  long long I0 = 1;
  long long init = 1;
  long long quotient;
  long long remainder;

  /* Compute the ideal decomposition, truncated to an integer, which
     minimizes the amount of communication. */
  rxs = (long long)sqrt((double)(Nx*Np)/(double)Ny);

  /* Constrain the decomposition */
  rx_max = Nx < Np ? Nx : Np; /* Require ry >= 1 and rx <= Nx */
  if(Ny < Np){ /* Require rx >= 1 and ry <= Ny */
    quotient = Np/Ny;
    remainder = Np%Ny;
    /* rx_min = the smallest integer >= Np/Ny */
    rx_min = quotient + (remainder > 0 ? 1 : 0);
  }
  else rx_min = 1;

  /* printf("rx_min = %d, rx_max = %d\n",rx_min, rx_max); */

  /* Constrain rxs to fall in this domain */
  rxs = rxs > rx_min ? rxs : rx_min;
  rxs = rxs < rx_max ? rxs : rx_max;

  /* Search down for a factor of Np */
  for(rx=rxs; rx>=rx_min; rx--){
    quotient = Np/rx;
    remainder = Np%rx;
    if(remainder == 0){
      rx0 = rx;
      ry0 = quotient;
      I0 = (rx0 - 1)*Ny + (ry0 - 1)*Nx;
      init = 0;
      break;
    }
  }

  /* Search up for a factor of Np */
  for(rx=rxs+1; rx<=rx_max; rx++){
    quotient = Np/rx;
    remainder = Np%rx;
    if(remainder == 0){
      ry = quotient;
      I = (rx - 1)*Ny + (ry - 1)*Nx;

      if(init || I < I0){
	rx0 = rx;
	ry0 = ry;
	I0  = I;
	init = 0;
      }
      break;
    }
  }

  if(init) return 1; /* Error locating a solution */

  /* printf("Minimum messaging decomposition has: rx = %d, ry = %d, I = %d\n",
     rx0, ry0, I0); */

  *pNGx = rx0;
  *pNGy = ry0;

  return 0;
}


/* ========================================================================== */


long long dom_decomp_3d(const long long Nx, const long long Ny, const long long Nz,
		  const long long Np, long long *pNGx, long long *pNGy, long long *pNGz){

  long long rx_min, rx_max, rx, ry, rz, I;
  long long rx0 = 1;
  long long ry0 = 1;
  long long rz0 = 1;
  long long I0 = 1;
  long long init=1;
  long long quotient;
  long long remainder;
  long long err, t, Npt;

  /* Constrain the decomposition */
  rx_max = Nx < Np ? Nx : Np; /* Require ry >= 1, rz >= 1 and rx <= Nx */
  /* Compute a global minimum constraint on rx. */
  t = (Ny < Np ? Ny : Np)*(Nz < Np ? Nz : Np); /* t = Max(ry)*Max(rz) */
  if(t < Np){ /* Require rx >= 1, ry <= Ny and rz <= Nz */
    quotient = Np/t;
    remainder = Np%t;
    /* rx_min = the smallest integer >= Np/t */
    rx_min = quotient + (remainder > 0 ? 1 : 0);
  }
  else rx_min = 1;

  /* printf("rx_min = %d, rx_max = %d\n",rx_min, rx_max); */

  for(rx = rx_min; rx <= rx_max; rx++){
    quotient = Np/rx;
    remainder = Np%rx;
    if(remainder == 0){
      Npt = quotient; /* Np for transverse (y,z) decomposition */

      err = dom_decomp_2d(Ny, Nz, Npt, &ry, &rz);
      if(err == 0){
	/* Now compute the amount of messaging */
	I = (rx - 1)*Ny*Nz + (ry - 1)*Nx*Nz + (rz - 1)*Nx*Ny;

	if(I < 0) continue; /* Integer Overflow */

	if(init || I < I0){
	  rx0 = rx;
	  ry0 = ry;
	  rz0 = rz;
	  I0  = I;
	  init = 0;
	  /* printf("I(rx = %d, ry = %d, rz = %d) = %d\n",rx,ry,rz,I); */
	}
      }
    }
  }

  if(init) return 1; /* Error locating a solution */

  *pNGx = rx0;
  *pNGy = ry0;
  *pNGz = rz0;

  return 0;
}

}// end namespace
/* ========================================================================== */

#ifdef DECOMP_TESTBED

#include <stdio.h>

/* ========================================================================== */

namespace PAMGEN_NEVADA {

long long dom_decomp_2d_serial(const long long Nx, const long long Ny,
			 const long long Np, long long *pNGx, long long *pNGy){

  long long rx, rx_min, rx_max, ry, I;
  long long rx0, ry0, I0, init=1;
  long long quotient;
  long long remainder;

  /* Constrain the decomposition */
  rx_max = Nx < Np ? Nx : Np; /* Require ry >= 1 and rx <= Nx */
  if(Ny < Np){ /* Require rx >= 1 and ry <= Ny */
    quotient = Np/Ny;
    remainder = Np%Ny;
    /* rx_min = the smallest integer >= Np/Ny */
    rx_min = quotient + (remainder > 0 ? 1 : 0);
  }
  else rx_min = 1;

  /* printf("rx_min = %d, rx_max = %d\n",rx_min, rx_max); */

  for(rx = rx_min; rx <= rx_max; rx++){
    quotient = Np/rx;
    remainder = Np%rx;
    if(remainder == 0){
      ry = quotient;

      /* Now compute the amount of messaging */
      I = (rx - 1)*Ny + (ry - 1)*Nx;

      if(init || I < I0){
	rx0 = rx;
	ry0 = ry;
	I0  = I;
	init = 0;
	/* printf("I(rx = %d, ry = %d) = %d\n",rx,ry,I); */
      }
    }
  }

  if(init) return 1; /* Error locating a solution */

  /* printf("Minimum messaging decomposition has: rx = %d, ry = %d, I = %d\n",
     rx0, ry0, I0); */

  *pNGx = rx0;
  *pNGy = ry0;

  return 0;
}


/* ========================================================================== */


long long dom_decomp_3d_serial(const long long Nx, const long long Ny, const long long Nz,
			 const long long Np, long long *pNGx, long long *pNGy, long long *pNGz){

  long long rx, rx_min, rx_max, ry, ry_min, ry_max, rz, I;
  long long rx0, ry0, rz0, I0, init=1;
  long long t, Npt;
  long long quotient;
  long long remainder;

  /* Constrain the decomposition */
  rx_max = Nx < Np ? Nx : Np; /* Require ry >= 1, rz >= 1 and rx <= Nx */
  /* Compute a global minimum constraint on rx. */
  t = (Ny < Np ? Ny : Np)*(Nz < Np ? Nz : Np); /* t = Max(ry)*Max(rz) */
  if(t < Np){ /* Require rx >= 1, ry <= Ny and rz <= Nz */
    quotient = Np/t;
    remainder = Np%t;
    /* rx_min = the smallest integer >= Np/t */
    rx_min = quotient + (remainder > 0 ? 1 : 0);
  }
  else rx_min = 1;

  /* printf("rx_min = %d, rx_max = %d\n",rx_min, rx_max); */

  for(rx = rx_min; rx <= rx_max; rx++){
    quotient = Np/rx;
    remainder = Np%rx;
    if(remainder == 0){
      Npt = quotient; /* Np for transverse (y,z) decomposition */
      ry_max = Ny < Npt ? Ny : Npt; /* Require rz >= 1 and ry <= Ny */
      if(Nz < Npt){ /* Require ry >= 1 and rz <= Nz */
        quotient = Npt/Nz;
        remainder = Npt%Nz;
	/* ry_min = the smallest integer >= Npt/Nz */
	ry_min = quotient + (remainder > 0 ? 1 : 0);
      }
      else ry_min = 1;

      /* printf("rx = %d, ry_min = %d, ry_max = %d\n",rx, ry_min, ry_max); */

      if(ry_min > ry_max)
	continue; /* No solution exists which satisfies the constraints */

      for(ry = ry_min; ry <= ry_max; ry++){
        quotient = Npt/ry;
        remainder = Npt%ry;
	if(remainder == 0){
	  rz = quotient;

	  /* Now compute the amount of messaging */
	  I = (rx - 1)*Ny*Nz + (ry - 1)*Nx*Nz + (rz - 1)*Nx*Ny;

	  if(I < 0) continue; /* Integer Overflow */

	  if(init || I < I0){
	    rx0 = rx;
	    ry0 = ry;
	    rz0 = rz;
	    I0  = I;
	    init = 0;
	    /* printf("I(rx = %d, ry = %d, rz = %d) = %d\n",rx,ry,rz,I); */
	  }
	}
      }
    }
  }

  if(init) return 1; /* Error locating a solution */

  *pNGx = rx0;
  *pNGy = ry0;
  *pNGz = rz0;

  return 0;
}


/* ========================================================================== */


long long main(void){

  long long Nx=512, Ny=1024, Nz = 256, Np = 8, NGx, NGy, NGz, Np_max;
  long long rx,ry,rz,I;
  long long err1, err2;

#if 0
  for(Np = 1; Np <= 2048; Np++){
    dom_decomp_2d(Nx, Ny, Np, &NGx, &NGy);

    printf("Nx = %d, Ny = %d, Np = %d, NGx = %d, NGy = %d\n",
	   Nx,Ny,Np,NGx,NGy);

    printf("Grids measure: Nx/NGx = nx = %e, Ny/NGy = ny = %e\n\n",
	   (double)Nx/(double)NGx, (double)Ny/(double)NGy);
  }
#endif

#if 1
  for(Nx = 1; Nx <= 48; Nx++){
    for(Ny = 1; Ny <= 48; Ny++){
      for(Np = 1; Np <= Nx*Ny; Np++){
	err1 = dom_decomp_2d(Nx, Ny, Np, &NGx, &NGy);
	err2 = dom_decomp_2d_serial(Nx, Ny, Np, &rx, &ry);

	if(err1 == 0 && err2 == 0){
	  if(rx != NGx || ry != NGy){
	    printf("Nx = %d, Ny = %d, Np = %d\n",Nx,Ny,Np);
	    printf("[std alg]: NGx = %d, NGy = %d\n",NGx,NGy);
	    printf("[linear alg]: NGx = %d, NGy = %d\n",rx,ry);
	  }
	  else if(NGx > Nx || NGy > Ny){
	    printf("Error: Nx = %d, Ny = %d, Np = %d, NGx = %d, NGy = %d\n",
		   Nx,Ny,Np,NGx,NGy);
	  }
	  /* else
	     printf("Nx=%d, Ny=%d, Np=%d, NGx=%d, NGy=%d\n",
	     Nx,Ny,Np,NGx,NGy);
	  */
	}
	else if(err1 + err2 == 1)
	  printf("err1 = %d, err2 = %d\n",err1, err2);
      }
    }
  }
#endif


#if 0
  for(Nx = 1; Nx <= 48; Nx++){
    for(Ny = 1; Ny <= 48; Ny++){
      for(Nz = 1; Nz <= 48; Nz++){
	Np_max = Nx*Ny*Nz;
	Np_max = Np_max < 100000 ? Np_max : 100000;
	for(Np = 1; Np <= Np_max; Np++){
	  err1 = dom_decomp_3d(Nx, Ny, Nz, Np, &NGx, &NGy, &NGz);
	  err2 = dom_decomp_3d_serial(Nx, Ny, Nz, Np, &rx, &ry, &rz);

	  if(err1 == 0 && err2 == 0){
	    if(rx != NGx || ry != NGy || rz != NGz){
	      printf("Nx = %d, Ny = %d, Nz = %d, Np = %d\n",Nx,Ny,Nz,Np);
	      printf("[std alg]: NGx = %d, NGy = %d, NGz = %d\n",NGx,NGy,NGz);
	      printf("[linear alg]: NGx = %d, NGy = %d, NGz = %d\n",rx,ry,rz);

	      I = (NGx - 1)*Ny*Nz + (NGy - 1)*Nx*Nz + (NGz - 1)*Nx*Ny;
	      printf("[std alg]: I = %d\n",I);

	      I = (rx - 1)*Ny*Nz + (ry - 1)*Nx*Nz + (rz - 1)*Nx*Ny;
	      printf("[linear alg]: I = %d\n",I);
	    }
	    else if(NGx > Nx || NGy > Ny || NGz > Nz){
	      printf("Error: Nx = %d, Ny = %d, Np = %d, NGx = %d, NGy = %d\n",
		     Nx,Ny,Np,NGx,NGy);
	    }
	    /* else
	       printf("Nx=%d, Ny=%d, Nz=%d, Np=%d, NGx=%d, NGy=%d, NGz=%d\n",
	       Nx,Ny,Np,NGx,NGy);
	    */
	  }
	  else if(err1 + err2 == 1)
	    printf("err1 = %d, err2 = %d\n",err1, err2);
	}
      }
    }
  }
#endif

  return 0;
}
}// end namespace  PAMGEN_NEVADA

/* ========================================================================== */

#endif /* DECOMP_TESTBED */

/* ========================================================================== */
