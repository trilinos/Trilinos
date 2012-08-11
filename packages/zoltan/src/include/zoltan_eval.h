/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef __ZOLTAN_EVAL_H
#define __ZOLTAN_EVAL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

struct Zoltan_Struct;


#define EVAL_LOCAL_SUM  0  /* calculated for nobj, obj_wgt, xtra_obj_wgt only */
#define EVAL_GLOBAL_SUM 1
#define EVAL_GLOBAL_MIN 2
#define EVAL_GLOBAL_MAX 3
#define EVAL_GLOBAL_AVG 4
#define EVAL_SIZE       5   /* must be last definition */

#define EVAL_MAX_XTRA_VWGTS 4
#define EVAL_MAX_XTRA_EWGTS 4

struct _eval_hg_struct{
  float obj_imbalance;          /* vertex number imbalance */
  float imbalance;              /* vertex weight imbalance */
  float cutl[EVAL_SIZE];        /* ConCut measure */
  float cutn[EVAL_SIZE];        /* NetCut measure */
  float nobj[EVAL_SIZE];        /* number of partition vertices */
  float obj_wgt[EVAL_SIZE];     /* partition vertex weights */
  float xtra_imbalance[EVAL_MAX_XTRA_VWGTS];
  float xtra_obj_wgt[EVAL_MAX_XTRA_VWGTS][EVAL_SIZE];
};

typedef struct _eval_hg_struct ZOLTAN_HG_EVAL;

struct _eval_graph_struct{
  float cuts[EVAL_SIZE];        /* The number of cut edges */
  float cut_wgt[EVAL_SIZE]  ;   /* The sum of the weights of the cut edges */
  float nnborparts[EVAL_SIZE];  /* The number of neighboring partitions */

  float obj_imbalance;          /* vertex number imbalance */
  float imbalance;              /* vertex weight imbalance */
  float nobj[EVAL_SIZE];        /* number of partition vertices */
  float obj_wgt[EVAL_SIZE];     /* partition vertex weights */
  float num_boundary[EVAL_SIZE];/* the number of objects with a remote neighbor */
              
  float xtra_imbalance[EVAL_MAX_XTRA_VWGTS];
  float xtra_obj_wgt[EVAL_MAX_XTRA_VWGTS][EVAL_SIZE];

  float xtra_cut_wgt[EVAL_MAX_XTRA_EWGTS][EVAL_SIZE];
};

typedef struct _eval_graph_struct ZOLTAN_GRAPH_EVAL;

struct _eval_balance_struct{
  float obj_imbalance;          /* vertex number imbalance */
  float imbalance;              /* vertex weight imbalance */
  float nobj[EVAL_SIZE];        /* number of partition vertices */
  float obj_wgt[EVAL_SIZE];     /* partition vertex weights */
              
  float xtra_imbalance[EVAL_MAX_XTRA_VWGTS];
  float xtra_obj_wgt[EVAL_MAX_XTRA_VWGTS][EVAL_SIZE];
};

typedef struct _eval_balance_struct ZOLTAN_BALANCE_EVAL;

int Zoltan_LB_Eval_Balance(struct Zoltan_Struct *zz, int print_stats, ZOLTAN_BALANCE_EVAL *eval);

int Zoltan_LB_Eval_Graph(struct Zoltan_Struct  *zz, int print_stats, ZOLTAN_GRAPH_EVAL *graph);

int Zoltan_LB_Eval_HG(struct Zoltan_Struct  *zz, int print_stats, ZOLTAN_HG_EVAL *hg);

int Zoltan_LB_Eval(struct Zoltan_Struct  *zz, int print_stats, 
                    ZOLTAN_BALANCE_EVAL *obj, ZOLTAN_GRAPH_EVAL *graph, ZOLTAN_HG_EVAL *hg);

void Zoltan_LB_Eval_Print_Graph(ZOLTAN_GRAPH_EVAL *graph);

void Zoltan_LB_Eval_Print_HG(ZOLTAN_HG_EVAL *hg);

void Zoltan_LB_Eval_Print_Balance(ZOLTAN_BALANCE_EVAL *lb);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
