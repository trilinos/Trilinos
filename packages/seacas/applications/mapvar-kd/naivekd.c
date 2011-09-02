/*
 * Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
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
/**************************************************************************/
/* Program to perform orthogonal range searches in a very naive k-d tree. */
/* In this implementation the nodes on any given level of the tree all    */
/* have the same discriminating dimension and the discriminator is chosen */
/* as NextDisc(i) = i+1 mod k.                                            */
/*                                                                        */
/* Later refinements employ much more sophisticated methods for the       */
/* selection of the discriminator.                                        */
/*                                                                        */
/* References:  J.L. Bentley "Multidimensional Binary Search Trees used   */
/* for Associative Searching.  ACM Sept. 1975 Vol. 18 No. 9.              */
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "naivekd.h"
/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP)  ((PP) = (TreePtr)malloc(sizeof(TreeRec)))
#define NEWTREE(PP)  if (TESTTREE(PP)==NULL) \
                         {printf("memory error\n");return;}


static TreePtr Root = NULL;

/*contains a pointer to an integer array of indices of points that were
  found in the rectangular query */

/***************************************************************************/
/* Determines if the treenode P falls inside the rectangular query         */
/* RectQuery.  If so, adds the array index of the point to the found       */
/* array.                                                                  */
/***************************************************************************/
void InRegion(TreePtr P,
	      int Dimension,
	      float *Points, int N, 
	      float *xmin, float *xmax,
	      int *found, int *count)
{
  int index, dc;

  index = P->Index;

  for (dc=0; dc < Dimension; dc ++) {
    if (Points[N*dc+index] < xmin[dc] || 
        Points[N*dc+index] > xmax[dc]) {
         return;
     }
  }
  /* P is in the region */
  found[(*count)++] = index+1;
}

/***************************************************************************/
/* Returns true iff the hyper-rectangle defined by bounds array B          */
/* intersects the rectangular query RectQuery.                             */
/***************************************************************************/
int BoundsIntersectRegion(float *B,
			  float *xmin, float *xmax,
			  int Dimension)
{
  int dc;

  for (dc=0; dc < Dimension; dc++) {
    if (B[2*dc] > xmax[dc] || B[2*dc+1] < xmin[dc]) {
      return(0);
    }
  }
  return(1);
}


/***************************************************************************/

void RangeSearch(TreePtr P,
		 float *Points, int N,
		 int Dimension,
		 float *xmin, float *xmax,
		 float *B,
		 int *found, int *count)
{
  int dc, disc;
  float *BHigh,*BLow;
  

  if (P==NULL) {printf("somehow a null pointer got sent here\n");}
  disc=P->Discrim;

  /* Check to see if the P is in the region */
  InRegion(P,Dimension,Points,N,xmin, xmax, found, count);

  BLow =  (float *)(malloc(2*Dimension*sizeof(float)));
  BHigh = (float *)(malloc(2*Dimension*sizeof(float)));

  if (BLow == NULL || BHigh == NULL) {
   printf("we have a memory error\n");
 }

  /* copy the region B into BLow, BHigh */
  for (dc=0; dc < 2*Dimension; dc++) {
    BLow[dc]  = B[dc];
    BHigh[dc] = B[dc];
  }

  /* Improve the Bounds for the subtrees */
  BLow[2*disc+1] = Points[N*disc+P->Index];
  BHigh[2*disc] = Points[N*disc+P->Index];

  if (P->Left != NULL && BoundsIntersectRegion(BLow, xmin, xmax, Dimension)) {
    RangeSearch(P->Left,Points,N,Dimension,xmin,xmax,BLow,found,count);
  }
  free(BLow);
  if (P->Right != NULL && BoundsIntersectRegion(BHigh, xmin, xmax, Dimension)) {
    RangeSearch(P->Right,Points,N,Dimension,xmin,xmax,BHigh,found,count);
  }
  free(BHigh);
}


/***************************************************************************/

void kdrectquery_(float *Points,
		  int *N,
		  int *Dimension,
		  float *xmin,
		  float *xmax,
		  int *found,
		  int *count)
{
float *B;
int dc;

  B =  (float *)(malloc(2* *Dimension*sizeof(float)));
  if (B == NULL) {
    printf("We have a memory problem\n");
  }

  for (dc =0; dc < *Dimension; dc++) {
    B[2*dc]   = xmin[dc];
    B[2*dc+1] = xmax[dc];
  }

  *count=0;
  RangeSearch(Root,Points,*N,*Dimension,xmin,xmax,B,found,count);
  /*  fprintf(stderr, "Found %d points\n", *count); */
  free(B);
}


/**************************************************************************/
int Successor(TreePtr P,
	      TreePtr Q,
	      int Dimension,
	      float *Points,
	      int N)
{
  int dc,disc;
  disc = Q -> Discrim;
  for (dc=0; dc < Dimension; dc++) {
    if (Points[N*((disc+dc) % Dimension)+P->Index] <
        Points[N*((disc+dc) % Dimension)+Q->Index]) 
      {
        return(-1);
      }
    if (Points[N*((disc+dc) % Dimension)+P->Index] >
        Points[N*((disc+dc) % Dimension)+Q->Index]) 
      {
        return(1);
      }
  }
  return(0);  /* they must be equal */

}


/***************************************************************************/
/* Inserts a node into the naive k-d tree */
void InsertNode(TreePtr P,
		float *Points,
		int N,
		int Dimension)
{
  int succ;
  TreePtr Son,Q;
  int index = P->Index;

  /* check for root */
  if (Root == NULL) {
    Root = P;
    P->Discrim = 0;
    P->Left=NULL;
    P->Right=NULL;
    return;
  }

  Son=Root;

  do {
    Q=Son;
    succ = Successor(P,Q,Dimension,Points,N);
    switch(succ) {
    case -1: Son=Q->Left;
      break;
    case  1: Son=Q->Right;
      break;
    case  0: return;  /* don't insert the point */
    }
  } while (Son != NULL);

  /* ASSERT: Q points to the leaf of the tree that needs to be added */
  if (succ==-1) {
    Q->Left = P;
  }
  else {
    Q->Right =P;
  }

  P->Discrim = (Q->Discrim + 1) % Dimension;
  P->Left = NULL;
  P->Right=NULL;  
}

/**************************************************************************/
void kdbuildtree_(float *Points,
		  int *NumPoints,
		  int *Dimension)
{
  int k;
  TreePtr Node;

  for (k=0; k < *NumPoints; k++) {
    NEWTREE(Node);
    Node->Index = k;
    InsertNode(Node,Points,*NumPoints, *Dimension);
  }
}

/***************************************************************************/
/*  Kills a kd-tree to avoid memory holes.   */
void KillTree(TreePtr P)
{

  if (P==NULL) {
    return;
  } /* just to be sure */
  if (P->Left != NULL) {
    KillTree(P->Left);
  }

  if (P->Right != NULL) {
    KillTree(P->Right);
  }

  free(P);
  P = NULL;
}

void kdkilltree_()
{
  KillTree(Root);
  Root = NULL;
}

