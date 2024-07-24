// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//Simple Util to help with test
#ifndef SHYLUBASKER_TEST_UTIL_HPP
#define SHYLUBASKER_TEST_UTIL_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

double myTime()
{
  struct timeval t1;
  double etime;
  double us = 1.0e6;
  
  gettimeofday(&t1,NULL);
  etime = t1.tv_sec*us + t1.tv_usec;
  return etime;
}

double totalTime(double start, double end)
{
  return (end-start)/1.0e6;
}

template <class Int, class Entry>
void readMatrix
(
 std::string fname, 
 Int &m, 
 Int &n, 
 Int &nnz,
 Int **col_ptr, 
 Int **row_idx,
 Entry **val
)
{
  std::ifstream inp_str;
  inp_str.open(fname, std::ios::in);
  
  std::string s;
  Int m_, n_, nnz_;
  m_ = n_ = nnz_ = 0;

  if(inp_str.is_open())
  {
    getline(inp_str, s);
    while ( inp_str.peek() == '%')
    { getline(inp_str, s); }

    inp_str >> m_; m = m_;
    inp_str >> n_; n = n_;
    inp_str >> nnz_; nnz = nnz_;

    //Check if need first alloc
    if(((*col_ptr) != NULL) && ((*col_ptr)[n_+1] < nnz_))
    {
      delete [] (*col_ptr); (*col_ptr) = NULL;
      delete [] (*row_idx); (*row_idx) = NULL;
      delete [] (*val); (*val) = NULL;
    }

    if((*col_ptr) == NULL)
    {
      *col_ptr = new Int[n_+1]();
      *row_idx = new Int[nnz_]();
      *val     = new Entry[nnz_]();
    }

    Int innz = 0;
    Int i_, j_;
    Entry v_;
    while(nnz_ > 0)
    {
      inp_str >> i_;
      (*row_idx)[innz] = i_-1;
      inp_str >> j_;
      (*col_ptr)[j_] = (*col_ptr)[j_]+1;
      inp_str >> v_;
      (*val)[innz] = v_;

      innz++;
      nnz_--;
    }//over nnz
    inp_str.close();
  }//is open

  for(Int k = 1; k < n+1; k++)
  {
    (*col_ptr)[k] = (*col_ptr)[k] + (*col_ptr)[k-1];
  }
}

template <class Int, class Entry>
void readVector(std::string fname, Int &n, Entry **x)
{
  std::string s;
  Int n_ = 0;
  Int nv_= 0;
  Entry v_ = (Entry) 0.0;
  std::ifstream inp_str;
  inp_str.open(fname, std::ios::in);
  
  if(inp_str.is_open())
  {
    getline(inp_str,s);
    while (inp_str.peek() == '%')
    { getline(inp_str,s); }

    inp_str >> n_;
    inp_str >> nv_; 
    if(nv_ > 1)
    { throw "TOO MANY VECTORS"; }

    //Check if alloc
    if(((*x)!= NULL) && (n < n_))
    {
      delete [] (*x); (*x) = NULL;
    }
    if((*x) == NULL)
    {
      *x = new Entry[n_]();
    }
    n = n_;
    Int ni = 0;
    while(n_ > 0)
    {
      inp_str >> v_;
      (*x)[ni] = (Entry) v_;
      //std::cout << "vect val" << (*x)[ni] << std::endl;
      n_--;
      ni++;
    }
    inp_str.close();
  }//if open
}//end readVector

template <class Int, class Entry>
Entry norm2(Int n, Entry x[])
{
  double sum = 0;
  for(Int i = 0; i < n; i++)
  {
    sum += x[i]*x[i];
  }
  sum = std::sqrt(sum);
  return sum;
}//norm2

template <class Int, class Entry>
void multiply
(
 Int m, 
 Int n, 
 Int col_ptr[],
 Int row_idx[], 
 Entry val[], 
 Entry x[], 
 Entry y[]
)
{
  for(Int i = 0; i < m; i++)
    {y[i] = (Entry) 0.0;}
  
  for(Int k = 0; k < n; k++)
  {
    for(Int i = col_ptr[k]; i < col_ptr[k+1]; i++)
    {
      const Int j = row_idx[i];
      y[j] += val[i]*x[k];
    }//over row
  }//over column
}//end multiply

template <class Int, class Entry>
void multiply_tr
(
 Int m, 
 Int n, 
 Int col_ptr[],
 Int row_idx[], 
 Entry val[], 
 Entry x[], 
 Entry y[]
)
{
  for(Int i = 0; i < m; i++)
    {y[i] = (Entry) 0.0;}
  
  // treat k as row in crs transpose, i.e. col in ccs non-tr
  for(Int k = 0; k < n; k++)
  {
    for(Int i = col_ptr[k]; i < col_ptr[k+1]; i++)
    {
      const Int j = row_idx[i]; // j is colid in crs transpose, i.e. rowid in ccs non-tr
      y[k] += val[i]*x[j];
    }//over row
  }//over column
}//end multiply

#endif
