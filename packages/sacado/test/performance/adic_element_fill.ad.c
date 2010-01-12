/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 01/12/10 10:38:06 by the version of   */
/*   ADIC compiled on  08/30/00 16:47:46                              */
/*                                                                    */
/*   ADIC was prepared as an account of work sponsored by an          */
/*   agency of the United States Government and the University of     */
/*   Chicago.  NEITHER THE AUTHOR(S), THE UNITED STATES GOVERNMENT    */
/*   NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, INCLUDING */
/*   ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS  */
/*   OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR */
/*   THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR  */
/*   PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE */
/*   PRIVATELY OWNED RIGHTS.                                          */
/*                                                                    */
/**********************************************************************/
#include "ad_deriv.h"
#include <math.h>
#include "adintrinsics.h"
#include <stdlib.h>
typedef struct  {
int  nqp;
int  nnode;
InactiveDouble  *w, *jac, **phi, **dphi;
int  *gid;
}
ElemData;
void   adic_element_fill(ElemData  *e,unsigned int  neqn,const DERIV_TYPE  *x,DERIV_TYPE  *u,DERIV_TYPE  *du,DERIV_TYPE  *f) {
unsigned int  var_0, var_1, var_2, var_3, var_4, var_5, var_6, var_7, var_8;
DERIV_TYPE  var_9;
double  adji_0;
    double  loc_0;
    double  loc_1;
    double  loc_2;
    double  loc_3;
    double  loc_4;
    double  loc_5;
    double  loc_6;
    double  loc_7;
    double  loc_8;
    double  loc_9;
    double  loc_10;
    double  loc_11;
    double  adj_0;
    double  adj_1;
    double  adj_2;
    double  adj_3;
    double  adj_4;

        static int g_filenum = 0;
        if (g_filenum == 0) {
            adintr_ehsfid(&g_filenum, __FILE__, "adic_element_fill");
        }
            for (unsigned int  qp = 0;     qp < e->nqp;     )    {
        for (unsigned int  eqn = 0;         eqn < neqn;         )        {
            {
                ad_grad_axpy_0(&(u[qp * neqn + eqn]));
                DERIV_val(u[qp * neqn + eqn]) = 0.0;
            }
            {
                ad_grad_axpy_0(&(du[qp * neqn + eqn]));
                DERIV_val(du[qp * neqn + eqn]) = 0.0;
            }
            for (unsigned int  node = 0;             node < e->nnode;             )            {
                {
                    loc_0 = DERIV_val(x[node * neqn + eqn]) * e->phi[qp][node];
                    loc_1 = DERIV_val(u[qp * neqn + eqn]) + loc_0;
                    ad_grad_axpy_2(&(u[qp * neqn + eqn]), 1.000000000000000e+00, &(u[qp * neqn + eqn]), e->phi[qp][node], &(x[node * neqn + eqn]));
                    DERIV_val(u[qp * neqn + eqn]) = loc_1;
                }
                {
                    loc_0 = DERIV_val(x[node * neqn + eqn]) * e->dphi[qp][node];
                    loc_1 = DERIV_val(du[qp * neqn + eqn]) + loc_0;
                    ad_grad_axpy_2(&(du[qp * neqn + eqn]), 1.000000000000000e+00, &(du[qp * neqn + eqn]), e->dphi[qp][node], &(x[node * neqn + eqn]));
                    DERIV_val(du[qp * neqn + eqn]) = loc_1;
                }
                var_2 = node++;
            }
            var_1 = eqn++;
        }
        var_0 = qp++;
    }
DERIV_TYPE  *s = malloc((e->nqp * neqn) *  sizeof (DERIV_TYPE ));
    for (unsigned int  qp = 0;     qp < e->nqp;     )    {
        for (unsigned int  eqn = 0;         eqn < neqn;         )        {
            {
                ad_grad_axpy_0(&(s[qp * neqn + eqn]));
                DERIV_val(s[qp * neqn + eqn]) = 0.0;
            }
            for (unsigned int  j = 0;             j < neqn;             )            {
                if (j != eqn)                 {
                    {
                        loc_0 = DERIV_val(s[qp * neqn + eqn]) + DERIV_val(u[qp * neqn + j]);
                        ad_grad_axpy_2(&(s[qp * neqn + eqn]), 1.000000000000000e+00, &(s[qp * neqn + eqn]), 1.000000000000000e+00, &(u[qp * neqn + j]));
                        DERIV_val(s[qp * neqn + eqn]) = loc_0;
                    }
                }
                var_5 = j++;
            }
            var_4 = eqn++;
        }
        var_3 = qp++;
    }
    for (unsigned int  node = 0;     node < e->nnode;     )    {
        for (unsigned int  eqn = 0;         eqn < neqn;         )        {
unsigned int  row = node * neqn + eqn;
            {
                ad_grad_axpy_0(&(f[row]));
                DERIV_val(f[row]) = 0.0;
            }
            for (unsigned int  qp = 0;             qp < e->nqp;             )            {
     DERIV_val(var_9) = exp(( DERIV_val(u[qp * neqn + eqn])));  /*exp*/
      adji_0 = DERIV_val(var_9);
                {
                    ad_grad_axpy_1(&(var_9), adji_0, &(u[qp * neqn + eqn]));
                }
                {
                    loc_0 = e->w[qp] * e->jac[qp];
                    loc_1 = e->jac[qp] * e->jac[qp];
                    loc_2 = 1.0 / loc_1;
                    loc_3 =  -loc_2;
                    loc_4 = loc_3 * DERIV_val(du[qp * neqn + eqn]);
                    loc_5 = loc_4 * e->dphi[qp][node];
                    loc_6 = DERIV_val(u[qp * neqn + eqn]) * DERIV_val(s[qp * neqn + eqn]);
                    loc_7 = DERIV_val(var_9) + loc_6;
                    loc_8 = e->phi[qp][node] * loc_7;
                    loc_9 = loc_5 + loc_8;
                    loc_10 = loc_0 * loc_9;
                    loc_11 = DERIV_val(f[row]) + loc_10;
                    adj_0 = e->phi[qp][node] * loc_0;
                    adj_1 = DERIV_val(u[qp * neqn + eqn]) * adj_0;
                    adj_2 = DERIV_val(s[qp * neqn + eqn]) * adj_0;
                    adj_3 = e->dphi[qp][node] * loc_0;
                    adj_4 = loc_3 * adj_3;
                    ad_grad_axpy_5(&(f[row]), 1.000000000000000e+00, &(f[row]), adj_4, &(du[qp * neqn + eqn]), adj_0, &(var_9), adj_2, &(u[qp * neqn + eqn]), adj_1, &(s[qp * neqn + eqn]));
                    DERIV_val(f[row]) = loc_11;
                }
                var_8 = qp++;
            }
            var_7 = eqn++;
        }
        var_6 = node++;
    }
    free(s);
}
void   AD_Init(int  arg0) {
    ad_AD_GradInit(arg0);

}
void   AD_Final() {
    ad_AD_GradFinal();

}
