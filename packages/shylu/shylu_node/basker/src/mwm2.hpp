// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MWM2_HPP
#define MWM2_HPP

//Maximum Weighted Matching
//This is based off Duff's mc64 
//-----------------------------------//

//Function list
// int mwm_init()       ... init fucntion
// int mwm()            ... main call function
// int mwm_trans        ... translates A-> \tilde{A}
// int mwm_carpaneto_row ... Carpaneto and Toth (1980) row init
// int mwm_carpaneto_col ... Carpaneto and Toth (1980) col init
// int mwm_diag_prod     ... Duff and Koster (1999) - mc64 Product of diag

//Struct 
// MWM_Q
// MWM_Q(size)
// MWM_Q.insert
// MQM_Q.empty
// MQM_Q.getMax()
//Queue structure based on Robert Sedgewick Algorith in C++ outline,
//Modified to to use top faster (Derigs and Metz (1986))

//Comeback and change
#include <float.h>
//#define INF      2000
#define INF      DBL_MAX

namespace mwm_order
{

  //-------------------------Data Structures-------------------//

  //Note:  In the future we should have a basker_q and use that 
  // with basker types
  // Still need to add Derigs and Metz adaptions
  //heap needs both an item and a key
  template <class Int, class Entry>
  struct MQM_Q
  {
    Entry *pq;
    Int   N,n;

    /*
    MWM_Q(Int size)
    {
      N = size;
      pq = new Entry[N+1];
      n = 0;
    }//end MWM_Q
    */

    //might need to add clear to the heap
    
    bool empty() const
    {
      return (N==0);
    }//end empty()

    void insert(Entry a)
    {
      pq[++N] = a;
      fix_up(pq, N);
    }//end insert
    
    Entry get_max()
    {
      exch(pq[1],pq[N]);
      fix_down(pq, 1, N-1);
      return pq[N--];
    }//end get_max
    
    void fix_up(Entry* a, Int k)
    {
      while(k > 1 &&  a[k/2] < a[k])
      {
        exch(a[k],a[k/2]);
        k = k/2;
      }
    }//end fix_up

    void fix_down(Entry *a, Int k, Int lN)
    {
      while(2*k <= lN)
      {
        Int j = 2*k;
        if(j < lN && a[j] < a[j+1])
        {
          j++;
        }
        if(!(a[k] < a[j]))
        {
          break;
        }
        exch(a[k], a[j]);
        k = j;	  
      }
    }//end fix_down

    void exch(Entry &a, Entry &b)
    {
      Entry temp = a;
      a = b;
      b = temp;
    }//end exch

  };// end MWM_Q
  

  template<class Int, class Entry>
  void mwm_heap_down2
  (
   Int i, 
   Int n,
   Int *Q,
   Entry *d,
   Int *L
  )
  {
    Int idum, pos, posk, qk;
    Entry di;
    
    const Int k = 2;
    
    pos = L[i];

    if(pos > 0)
    {
      di = d[i];

      for(idum = 0; idum < n; idum++)
      {
        posk = pos/k;
        qk   = Q[posk];

        if(di <= d[qk])
        {
          break;
        }

        Q[pos] = qk;
        L[qk]  = pos;
        pos    = posk;

        if(pos-1 < 0)
        {
          break;
        }

      }//for -idum
    }//if - pos > 0
    
    Q[pos] = i;
    L[i]   = pos;
  }//end mwm_hep_down_2


  template <class Int, class Entry>
  void mwm_heap_down
  (
   Int i,
   Int n, 
   Int *Q,
   Entry *D,
   Int *L
  )
  {
    //Note: We could really clean up the logic here

    Int k = 2;    
    Int pos = L[i];
    if(pos != 0)
    {
      Entry di = D[i];
      for(Int idum = 0; idum < n; idum++)
      {
        Int posk = pos/k;
        Int qk   = Q[posk];

        if(di >= D[qk])
        {
          break;
        }

        #ifdef MATCH_DEBUG
        printf("heap down pos: %d \n", pos);
        #endif

        Q[pos] = qk;
        L[qk]  = pos;
        pos    = posk;

        if(pos <=0)
        {
          break;
        }
      }//for-idum
    }//if-pos!=0
    
    #ifdef MATCH_DEBUG
    printf("heap down pos2: %d \n", pos);
    #endif

    Q[pos] = i;
    L[i]   = pos;
    
    //Done!
  }//end mwm_heap_down


  template <class Int , class Entry>
  void mwm_heap_del_root2
  (
   Int &qlen,
   Int n,
   Int *Q,
   Entry *d,
   Int *L
  )
  {
    Int i, idum, pos, posk;
    Entry dk, dr, di;
    Int k = 2;

    i  = Q[qlen-1];
    di = d[i];
    qlen = qlen-1;
    pos = 0;

    for(idum = 0; idum < n; idum++)
    {
      posk = k*pos;
      if(posk < qlen)
      {
        dk = d[Q[posk]];
        //-1
        if(posk < qlen)
        {
          dr = d[Q[posk+1]];

          if(dk < dr)
          {
            posk = posk+1;
            dk = dr;
          }

        }

        if(di >= dk)
        {
          break;
        }
      }
      else
      {
        break;
      }

      Q[pos] = Q[posk];
      L[Q[pos]] = pos;
      pos = posk;
    }
    
    Q[pos] = i;
    L[i] = pos;
  }//end mwm_head_del_root2

  
  template <class Int, class Entry>
  void mwm_heap_del_root
  (
   Int &qlen,
   Int n,
   Int *Q,
   Entry *D,
   Int *L
  )
  {
    const Int k = 2;

    Int i = Q[qlen];
    Entry di = D[i];
    
    qlen--;
    Int pos = 0;
    
    for(Int idum = 0; idum < n; idum++)
    {
      Int posk = k*pos;
      if(posk > qlen)
      {
        break;
      }
      Entry dk = D[Q[posk]];
      if(posk < qlen)
      {
        Entry dr = D[Q[posk+1]];
        if(dk > dr)
        {
          posk++;
          dk = dr;
        }
      }
      if(di <= dk)
      {
        break;
      }
      Q[pos]    = Q[posk];
      L[Q[pos]] = pos;
      pos       = posk; 	
    }//for--idum(1:N)
    
    Q[pos] = i;
    L[i] = pos;

    //Done
  }//end mwm_heap_del_root


  template <class Int, class Entry>
  void mwm_heap_extract_reheap2
  (
   Int pos0,
   Int &qlen,
   Int n,
   Int *Q,
   Entry *d,
   Int *L
  )
  {
    Int i, idum,pos, posk, qk;
    Entry dk, dr, di;
    
    const Int k = 2;
    
    if(qlen-1 == pos0)
    {
      qlen--;
      return;
    }

    i  = Q[qlen-1];
    di = d[i];
    qlen = qlen-1;
    pos  = pos0;
    
    if(pos > 0)
    {
      for(idum = 0; idum < n; idum++)
      {
        posk = pos/k;
        qk = Q[posk];

        if(di <= d[qk])
        {
          break;
        }

        Q[pos] = qk;
        L[qk]  = pos;
        pos = posk;

        if(pos <=0)
        {
          break;
        }
      }

      Q[pos] = i;
      L[i] = pos;

      if(pos != pos)
      {
        return;
      }

      for(idum = 0; idum < n; idum++)
      {
        posk = k*pos;

        if(posk >= qlen)
        {
          break;
        }

        dk = d[Q[posk]];

        if(posk < qlen-1)
        {
          dr = d[Q[posk+1]];
          if(dk < dr)
          {
            posk = posk+1;
            dk = dr;
          }
        }

        if(di >= dk)
        {
          break;
        }

        qk = Q[posk];
        Q[pos] = qk;
        L[qk] = pos;
        pos   = posk;
      }
    }
    Q[pos] = i;
    L[i]   = pos;
    //Done
  }//end mwm_heap_extract_reheap2


  template <class Int, class Entry>
  void mwm_heap_extract_reheap
  (
   Int pos0,
   Int &qlen,
   Int n,
   Int *Q,
   Entry *D,
   Int *L
  )
  {
    Int k = 0; // NDE: why is this zero and not 2, like above?
    
    if(qlen == pos0)
    {
      qlen--;
      return;
    }
    
    //? qlen-1?
    Int i    = Q[qlen];
    Entry di = D[i];
    qlen--;
    Int   pos = pos0;

    if(pos >= 0)
    {
      //go down
      for(Int idum = 0; idum < n; idum++)
      {
        Int posk = pos/k; // NDE: will this lead to div by zero? k init to 0 above and not updated...
        Int qk   = Q[posk];

        if(di >= D[qk])
        {
          break;
        }

        Q[pos] = qk;
        L[qk]  = pos;
        pos    = posk;

        if(pos <= 0)
        {
          break;
        }
      }

      Q[pos] = i;
      L[i] = pos;

      if(pos != pos0)
      {
        return;
      }

      //go up
      for(Int idum = 0; idum <n; idum++)
      {
        Int posk = k*pos; // NDE: will this be zero?

        if(posk > qlen)
        {
          break;
        }

        Entry dk = D[Q[posk]];

        if(posk < qlen)
        {
          Entry dr = D[Q[posk+1]];
          if(dk > dr) 
          {
            posk++;
            dk = dr;
          }
        }

        if(di <= dk)
        {
          break;
        }

        Int qk     = Q[posk];
        Q[pos] = qk; 
        L[qk]  = pos;
        pos    = posk;
      }//for-idum

    } //if pos>=0

    Q[pos] = i;
    L[i]   = pos;

    //Done
  }//end mwm_heap_extract_reheap


  //---------------------------Functions-----------------///
  //Function list
  // int mwm_init()       ... init fucntion
  // int mwm()            ... main call function
  // int mwm_trans        ... translates A-> \tilde{A}
  // int mwm_carpaneto_row ... Carpaneto and Toth (1980) row init
  // int mwm_carpaneto_col ... Carpaneto and Toth (1980) col init
  // int mwm_diag_prod     ... Duff and Koster (1999) - mc64 Product of diag

   
  template<class Int, class Entry>
  int mwm_bn_init
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry  *val,
   Int *pr, 
   Int *L,
   //Entry *d, 
   double *d, 
   Int *iperm, 
   Int *jperm, 
   Int &num, 
   double &bv
  )
  {
    Int i,ii,i0,j,jj, k;
    Int kk, kk1, kk2;
//    Entry a0, ai;
    double a0, ai;
    bv = (double) INF;

    //Init used values
    i0 = -1;
    num = 0;
    for(k = 0; k < n; k++)
    {
      iperm[k] = -1;
      jperm[k] = -1;
      pr[k]    = col_ptr[k];
      //d[k]     = (Entry) 0;
      d[k]     = (double) 0;
    }
    
    //Scan over column nodes
    for(j=0; j<n; j++)
    {
      //a0 = (Entry) -1.0;
      a0 = (double) -1.0;
      //For each column node, 
      for(k=col_ptr[j]; k<col_ptr[j+1]; k++)
      {
        i  = row_idx[k];
        ai = abs(val[k]);
        if(ai > d[i])
        {
          d[i] = ai;
        }
        if(jperm[j] != -1)
        {
          continue;
        }
        if(ai >= bv)
        {
          a0 = bv;
          if(iperm[i] != -1)
          {
            continue;
          }
          jperm[j] = i;
          iperm[i] = j;

          num++;
        }
        else
        {
          if(ai <= a0)
          {
            continue;
          }
          a0 = ai;
          i0 = i;
        }
      }//for-k, row nodes
      if((a0 != ((double)(-1.0))) && (a0 < bv))
      {
        bv = a0;
        if(iperm[i0] != -1)
        {
          continue;
        }
        iperm[i0] = j;
        jperm[j]  = i0;

        num++;
      }
    }//for-j

    //Update BV with smallest of all largest max |val| row
    for(i = 0; i < n; i++)
    {
      bv = min(bv,d[i]);
    }

    if(num == n)
    {
      return 0;
    }

    //if not shift 
    for(j =0; j < n; j++)
    {
      if(jperm[j]!=-1)
      {
        continue;
      }
      for(k=col_ptr[j]; k<col_ptr[j+1]; k++)
      {
        i  = row_idx[k];
        ai = abs(val[k]);
        bool d_break = false;

        if(ai < bv)
        {
          continue;
        }

        if(iperm[i] == -1)
        {
          num++;
          jperm[j] = i;
          iperm[i] = j;
          pr[j]    = k+1;
          break;
        }

        jj  = iperm[i];

        kk1 = pr[jj];
        kk2 = col_ptr[jj+1];

        if(kk1 >= kk2)
        {
          continue;
        }

        for(kk=kk1; kk<kk2; kk++)
        {
          ii = row_idx[kk];
          if(iperm[ii] != -1)
          {
            continue;
          }
          if(abs(val[kk])>=bv)
          {
            jperm[jj] = ii;
            iperm[ii] = jj;
            pr[jj]    = kk+1;
            num++;
            jperm[j] = i;
            iperm[i] = j;
            pr[j] = k+1;

            d_break = true;
            break;
          }//if- 
        }//for-kk

        if(d_break == true)
        {
          d_break = false;
          break;
        }

        pr[jj] = kk2 +1;
      }//for-k
    }//for-j
    //DONE
    return 0;
  }//end mwm_bn_init


  template<class Int, class Entry>
  int mwm_bn
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry  *val,
   Int *pr, 
   Int *L,
   //Entry *d, 
   double *d, 
   Int *iperm, 
   Int *jperm, 
   //Int &num, 
   //Entry &bv,
   Int &num, 
   double &bv
  )
  {
    Int i, i0;
    Int j, jj;
    Int k, kk;;
    Int jord, jdum, idum;
    Int qlen, low, up;
    Int q0;
    double dq0;
    //Entry dq0;

    double dnew, di;
    double csp;
    //Entry dnew, di;
    //Entry csp;
    Int isp, jsp;
    
    Int lpos;

    //Entry MINONE = (Entry) -1.0;
    double MINONE = (double) -1.0;

    Int *Q = new Int[n+1];
   
    //Debug
    #ifdef MWM_DEBUG
    printf("iperm:\n");
    for(Int jj = 0; jj < n; jj++)
    {
      printf("%d, ", iperm[jj]);
    }
    printf("\n");
    printf("jperm: \n");
    for(Int jj = 0; jj < n; jj++)
    {
      printf("%d, ", jperm[jj]);
    }
    printf("\n");
    #endif
    
    jsp = -1;
    isp = -1;
    //prep
    for(i = 0; i < n; i++)
    {
      d[i] = MINONE;
      L[i] = -1;
    }

    for(jord = 0; jord < n; jord++)
    {
      if(jperm[jord] != -1)
      {
        continue;
      }

      qlen = 0;
      low  = n;
      up   = n;

      csp   = MINONE;
      j     = jord;
      pr[j] = -1;


      for(k=col_ptr[j]; k < col_ptr[j+1]; k++)
      {
        i    = row_idx[k];
        dnew = abs(val[k]);

        if(csp >= dnew)
        {
          continue;
        }
        if(iperm[i] == -1)
        {
          csp = dnew;
          isp = i;
          jsp = j;

          if(csp >= bv)
          {
            goto L160;
          }//if-csp>=bv
        }
        else
        {
          d[i] = dnew;
          if(dnew >= bv)
          {
            low--;
            Q[low] = i;
          }
          else
          {
            L[i] = qlen;
            qlen++;
            mwm_heap_down2(i,n,Q,d,L);
          }
          jj = iperm[i];
          pr[jj] = j;
        }// end else
      }//for-k

      for(jdum = 0; jdum < num; jdum++)
      {
        bool d_break = false;

        if(low == up)
        {
          if(qlen == 0)
          {
            goto L160;
          }

          i = Q[0];

          if(csp >= d[i])
          {
            goto L160;
          }

          bv = d[i];

          for(idum = 0; idum < n; idum++)
          {
            mwm_heap_del_root2(qlen,n,Q,d,L);
            L[i] = -1;
            low--;
            Q[low] = i;

            if(qlen == 0)
            {
              break;
            }

            i = Q[0];

            if(d[i] != bv)
            {
              break;
            }
          }//for-idum
        }//if-low==up

        up--;
        q0   = Q[up];
        dq0  = d[q0];
        L[q0] = up;

        j = iperm[q0];

        for(k = col_ptr[j]; k < col_ptr[j+1]; k++)
        {
          i = row_idx[k];

          if(L[i] >= up)
          {
            continue;
          }

          dnew = min(dq0, abs(val[k]));
          
          if(csp >= dnew)
          {
            continue;
          }
          if(iperm[i] == -1)
          {
            csp = dnew;
            isp = i;
            jsp = j;
            if(csp >= bv)
            {
              goto L160;
            }
          }
          else
          {
            di = d[i];
            if((di >= bv) || (di>=dnew))
            {
              continue;
            }

            d[i] = dnew;

            if(dnew >= bv)
            {
              if(di != MINONE)
              {
                lpos = L[i];
                mwm_heap_extract_reheap2(lpos,
                    qlen, n, Q, d, L);

              }//if- di!= minone
              L[i] = -1;
              low--;
              Q[low] = i;
            }
            else
            {
              if(di == MINONE)
              {
                L[i] = qlen;
                qlen++;
              }

              mwm_heap_down2(i,n,Q,d,L);
            }//if-else

            jj = iperm[i];
            pr[jj] = j;
          }//elws			 
        }//for-k
        if(d_break == true)
        {
          break;
        }
      }//for-jdum (150)

L160:
      //(160)
      if(csp != MINONE)
      {
        bv = min(bv,csp);
        num++;
        i = isp;
        j = jsp;

        for(jdum = 0; jdum < num+1; jdum++)
        {
          i0 = jperm[j];
          jperm[j] = i;
          iperm[i] = j;
          //printf("Backtrack iperm[%d] = %d \n",
          //     i, iperm[i]);
          j = pr[j];
          if(j == -1)
          { 
            break;
          }
          i = i0;
        }
      }//if-csp != minone 

      for(kk=up; kk < n; kk++)
      {
        i = Q[kk];
        d[i] = MINONE;
        L[i] = -1;
      }
      //up or up-1
      for(kk=low; kk < up; kk++)
      {
        i = Q[kk];
        d[i] = MINONE;
      }

      for(kk = 0; kk < qlen; kk++)
      {
        i = Q[kk];
        d[i] = MINONE;
        L[i] = -1;
      }

    }//for jord

    if(num!=n)
    {
      //printf("Struct singluar\n");
      for(j = 0; j < n; j++)
      {
        jperm[j] = -1;
      }
      k = 0;
      for(i=0; i < n; i++)
      {
        if(iperm[i]==-1)
        {
          k++;
          pr[k] = i;
        }
        else
        {
          j = iperm[i];
          jperm[j] = i;
        }
      }
      for(i = 0; i < n; i++)
      {
        if(jperm[i] != -1)
        {
          continue;
        }
        k++;
        jdum = pr[k];
        iperm[jdum] = -i;
      }
    }//if(n != num)

    delete [] Q;
    return 0;
  }//end mwm_bn


  template <class Int, class Entry>
  int mwm_init()
  { return -1;}
  

  template <class Int, class Entry>
  int mwm
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry *val, 
   Int *perm,
   Int &num
  )
  {
    double *d    = new double[n];
    Int   *jperm = new Int[n];
    Int   *iperm = new Int[n];
    Int   *L     = new Int[n];
    Int   *pr    = new Int[n];

    //Entry bv = 0;
    double bv = 0;

    mwm_bn_init(n,nnz, 
		col_ptr, row_idx, val,
		pr, L, d,
		iperm, jperm, 
		num, bv);

    #ifdef MATCH_DEBUG
    printf("\n");
    printf("Bottleneck init done.  num: %d n: %d \n", num, n);
    printf("\n");
    #endif

    for(Int i = 0; i < n; i++)
    {
      perm[i] = iperm[i];
    }

    #ifdef MATCH_DEBUG
    FILE *fp;
    fp = fopen("bn_init.txt", "w");
    printf("BN init perm: \n");
    for(Int i = 0; i < n; i++)
    {
      fprintf(fp, "%ld \n", iperm[i]);
      //printf("%ld, ", iperm[i]);
    }
    fclose(fp);
    printf("\n");
    #endif

    if(num == n)
    {
      delete [] d;
      delete [] jperm;
      delete [] iperm;
      delete [] L;
      delete [] pr;
      return 0;
    }

    mwm_bn(n,nnz,
	   col_ptr, row_idx, val,
	   pr, L, d, 
	   iperm, jperm,
	   num, bv);
    
    #ifdef MATCH_DEBUG
    printf("\n");
    printf("Bottleneck done. num: %d \n", num);
    printf("\n");
    #endif

    for(Int i = 0; i < n; i++)
    {
      perm[i] = iperm[i];
    }

    #ifdef MATCH_DEBUG
    printf("BN perm: \n");
    for(Int i = 0; i < n; i++)
    {
      printf("%d, ", iperm[i]);
    }
    printf("\n");
    #endif

    delete [] d;
    delete [] jperm;
    delete [] iperm;
    delete [] L;
    delete [] pr;
    return 0;
  }//end mwm()


  //Main calling driver function
  template <class Int, class Entry>
  int mwm_prod
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry *val, 
   Int *perm, 
   Int &num
  )
  {
    #ifdef MATCH_DEBUG
    printf("\n");
    printf("---------------MWM CALLED-----------------");
    printf("\n");
    #endif

    Entry *U          = new Entry[n];
    Entry *d          = new Entry[n];
    Int   *jperm      = new Int[n];
    Int   *iperm      = new Int[n];
    Int   *L          = new Int[n];
    Int   *pr         = new Int[n];
    Entry *min_val    = new Entry[nnz+1];
    mwm_tras(n, nnz, col_ptr, row_idx, val, min_val);

    #ifdef MATCH_DEBUG
    printf("\n");
    printf("---------After val translation-------");
    printf("\n");
    printf("min_val: \n");
    for(Int i = 0; i < nnz; i++)
    {
      printf("%e, ", min_val[i]);
    }
    printf("\n");
    #endif

    //---------------Starting matching---///
    //Rows
    num =0;
    //printf("num: %d \n", num);

    mwm_carpaneto_row(n,nnz,col_ptr,row_idx,min_val,
		      pr,L,U,d,
		      iperm,jperm,num);
    
    #ifdef MATCH_DEBUG
    printf("\n");
    printf("------------Done row: %d -----------\n", num);
    printf("\n");
    printf("Dual Val U: \n");
    for(Int i = 0; i < n; i++)
    {
      printf("%e, ", U[i]);
    }
    printf("\n");
    #endif


    if(num == n)
    {
      delete [] U;
      delete [] d;
      delete [] jperm;
      delete [] iperm;
      delete [] L;
      delete [] pr;
      delete [] min_val;
      return 0;
    }
    
    //Columns
    mwm_carpaneto_col(n,nnz,col_ptr,row_idx,min_val,
		      pr, L, U, d, 
		      iperm, jperm,num);

    #ifdef MATCH_DEBUG
    printf("\n");
    printf("--------------Done col: %d -------------\n", num);
    printf("\n");
    printf("Dual Val U: \n");
    for(Int i = 0; i < n; i++)
    {
      printf("%e, ", U[i]);
    }
    printf("\n");
    #endif
 
    if(num == n)
    {
      delete [] U;
      delete [] d;
      delete [] jperm;
      delete [] iperm;
      delete [] L;
      delete [] pr;
      delete [] min_val;
      return 0;
    }

    #ifdef MATCH_DEBUG
    printf("\n MATCH PERM: \n");
    for(Int i=0; i < n; i++)
    {
      printf("%d \n", iperm[i]);
    }
    printf("\n");
    #endif

    //----------------Do Matching ------//
    //printf("calling diag prod\n");
    mwm_diag_prod(n, nnz, col_ptr, row_idx, min_val,
		  pr, L, U, d, 
		  iperm, jperm, num);

    //Until done with debug, copy perm
    //for(Int i = 0; i < n; i++)
    //{
      //perm[i] = iperm[i];
    //}

    #ifdef MATCH_DEBUG
    printf("\n MATCH PERM: \n");
    for(Int i=0; i < n; i++)
    {
      printf("%d \n", iperm[i]);
    }
    printf("\n");
    #endif
		  
    delete [] U;
    delete [] d;
    delete [] jperm;
    delete [] iperm;
    delete [] L;
    delete [] pr;
    delete [] min_val;
    return 0;
  }//end mwm()


  //Converts to do a mini opt problem
  template <class Int, class Entry>
  int mwm_tras
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry *val,  
   Entry *min_val
  )
  {
    for(Int k = 0; k < n; k++)
    {
      Entry fact = 0.0;
      //Find max value in abs in each column
      for(Int j = col_ptr[k]; j < col_ptr[k+1]; j++)
      {
        if(abs(val[j]) > fact)
        {
          fact = abs(val[j]);
        }
      }//end nnz in column

      for(Int j = col_ptr[k]; j < col_ptr[k+1]; j++)
      {
        //come back to the left shift op
        min_val[j]= fact - abs(val[j]);
      }
    }//end over all columns    

    return 0;
  }//end mwm_trans()


  //This finds an inital matching based off the dual value for row
  // $$U(i) =min_{j \in ROW(i)} c_{ij} \;  \forall i \in V_{r}
  template <class Int, class Entry>
  int mwm_carpaneto_row
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry  *val,
   Int *pr, 
   Int *L,
   Entry   *U, 
   Entry *d, 
   Int *iperm, 
   Int *jperm, 
   Int &num
  )
  {
    //printf("-----MWM  Row called -----\n");
    //Init values
    num = 0;
    for(Int k = 0; k < n; k++)
    {
      U[k]     = (Entry) INF;
      d[k]     = (Entry) 0;
      iperm[k] = -1;
      jperm[k] = -1;
      pr[k]    = col_ptr[k];
      L[k]     = -1;
    }//end for over all columns

    //Preload all U, though do not select yet
    for(Int k = 0; k < n; k++)
    {
      for(Int j = col_ptr[k]; j < col_ptr[k+1]; j++)
      {
        Int i = row_idx[j];
        if(val[k] < U[i])
        {
          U[i]      = val[k];
          iperm[i]  = k;
          L[i]      = j;
        }
      }//end over all nnz in column
    }//end over all columns
   
    //Now use the preloads to find a permuation
    for(Int i = 0; i < n; i++)
    {
      Int j = iperm[i];
      //if row permuted
      if(j != -1)
      {
        iperm[i] = 0;
        //If that column has not been match yet
        if(jperm[j] == -1)
        {
          //if not dense colum
          //we use the same standards are mc64
          if(((col_ptr[j+1]-col_ptr[j]) < (n/10)) ||
              (n < 50))
          {
            num+=1;
            iperm[i] = j;
            jperm[j] = L[i];
          }//not dense
        }//column not match
      }//row populated and not matched
    }//end over all rows
    
    return 0;
  }//end mwm_carpaneto_row()


  //This finds an inital matching based off the dual value for row
  // $$v_{j} = min_{i \in COL(j)} (c_{ij} - u_{i}) \forall j \in V_{c} 
  template <class Int, class Entry>
  int mwm_carpaneto_col
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry  *val,
   Int *pr, 
   Int *L,
   Entry   *U, 
   Entry *d, 
   Int *iperm, 
   Int *jperm, 
   Int &num
  )
  {
    //printf("-------MWM col: %d -----------\n", num); 

    Int k, i,ii,j,jj;
    Int k0,k1, k2;
    Int kk,kk0, kk1, kk2;
    Int i0;

    Entry di, vj;


    for(j = 0; j < n; j++)
    {
      //if column has not been matched
      bool found_replace = false;
      if(jperm[j] == -1)
      {
        k1   = col_ptr[j]; //should use const
        k2   = col_ptr[j+1];

        if(k1 > k2)
        {
          continue;
        }

        i0   = row_idx[k1];
        vj = val[k1];
        k0 = k1;

        for( k = k1+1; k < k2; k++)
        {
          i    = row_idx[k];
          di   = val[k] - U[i];

          //Might be a better way for 
          //this flow
          if(di > vj)
          {
            continue;
          }

          //if min than already
          if((di <= vj)||(di==INF))
          {
            if(di == vj)
            {
              //if i already mactch
              //or i0 is not matched
              if((iperm[i]  != -1) ||
                  (iperm[i0]  == -1))
              {
                //continue over all 
                //possible
                continue;
              }
            }
            vj = di;
            i0 = i;
            k0 = k;
          }//if should consider value
        }//over all possible nnz


        //vj should be the min.
        //i0 the min row_idx
        //
        //now assign if not already assin
        //if already assign try to find alt
        d[j] = vj;
        k = k0;
        i = i0;

        //Already assigned need to
        //try to find alt
        //bool found_replace = false;
        if(iperm[i] != -1)
        {
          for(k = k0; k <k2; k++)
          {
            i = row_idx[k];

            //if an alt value
            if((val[k]-U[i])<= vj)
            {
              jj  = iperm[i];

              kk1 = pr[jj];
              kk2 = col_ptr[jj+1];

              if(kk1 > kk2)
              {
                continue;
              }

              for(kk = kk1;
                  kk < kk2;
                  kk++)
              {
                ii = row_idx[kk];

                if(iperm[ii] != 0)
                {

                  if(val[kk]-U[ii] <= d[jj])
                  {
                    #ifdef MATCH_DEBUG
                    printf("MWM COL replace\n");
                    #endif

                    jperm[jj] = kk;
                    iperm[ii] = jj;
                    pr[jj] = kk+1;
                    found_replace = true;
                    break;
                  }//if right value
                }//if new is notmatched
              }//for --
              //over all new possible
              //maybe we can find bter
            }//if possible good value

            if(found_replace == true)
            {
              break;
            }
            pr[jj] = kk2+1;
          }//for --over all other column
        }//if found is already matched
        else
        {
          found_replace = true;
          #ifdef MATCH_DEBUG
          printf("MWM COL found not matched\n");
          #endif
        }
      }//not permuted
      else
      {
        continue;
      }

      //We can make this logic better
      if(found_replace == true)
      {
        #ifdef MATCH_DEBUG
        printf("MWM COL, add: i: %d j: %d \n", 
            i, j);
        #endif 

        //might have to make k for all
        //might have to make i for all
        num++;
        //printf("i: %d  j: %d \n", i, j);
        jperm[j] = k;
        iperm[i] = j;
        pr[j]    = k+1;
        found_replace = false;
      }
    }//for --over all columns
    
    return 0;
  }//end mwm_carpaneto_col()

  template <class Int, class Entry>
  int mwm_diag_prod
  (
   Int n, 
   Int nnz,
   Int *col_ptr, 
   Int *row_idx,
   Entry  *val,
   Int *pr, 
   Int *L,
   Entry *U, 
   Entry *d, 
   Int *iperm, 
   Int *jperm, 
   Int &num
  )
  {
    Int i, k;
    Int isp, jsp;
    Int jord;

    Entry dnew;

    Int *Q = new Int[n+1]; //Q
    Int *out = new Int[n+1];

    //reinit varaibles
    for(k = 0; k < n; k++)
    {
      d[k] = INF;
      L[k] = -1;
    }//end for over all rows

    //Each loop is similar to Dijkstra's alg
    //Solving single source shortest path problem
    for(jord = 0; jord < n; jord++)
    {
      if(jperm[jord] != -1)
      {
        continue;
      }

      Entry dmin = INF; //len of the shortest path in tree
      Entry csp  = INF; //cost of augmented path

      //Q stuff
      Int   qlen = 0;
      Int   low  = n;
      Int   up   = n;

      //if(jperm[jord]==-1)
      {
        Int j = jord;
        pr[j] = -1; // the root for j


        //Scan columns of J
        for( k = col_ptr[j]; k < col_ptr[j+1]; k++)
        {
          i = row_idx[k];
          dnew = val[k] - U[i];

          //for debuggin
          if(dnew < 0)
          {
            printf("\n dnew < 0 \n");
            printf("j: %d i: %d k: %d val: %e U: %e dnew: %e \n", j, i, k, val[k], U[i], dnew);

            delete [] Q;
            delete [] out;
            return -1;
          }

          if(dnew >= csp)
          {
            continue;
          }
          //if(dnew < csp)
          {
            //if row has not been matched
            if(iperm[i] == -1)
            {
              csp = dnew;
              isp = k;
              jsp = j;
            }
            else //if matched
            {
              if(dnew < dmin)
              {
                dmin = dnew;
              }
              d[i] = dnew;
              Q[qlen++] = k;
            }
          }//if new value is less
        }//for-(k) over nnz in J

        #ifdef MATCH_DEBUG
        printf("MWM, col: %d qlen: %d dmin:%e csp: %e \n",
            j, qlen, dmin, csp);
        #endif

        //heapify what we have found
        //Inital heap Q      rows held in Q(1:qlen)
        Int q0 = qlen;
        qlen = 0;
        for(Int kk = 0; kk < q0; kk++)
        {
          k  = Q[kk];
          i  = row_idx[k];
          //if distance is larger than current path
          //don't add
          if(csp <= d[i])
          {
            d[i] = INF;
            continue;
          }
          //Add this to the top Q
          if(d[i] <= dmin)
          {
            low--;
            Q[low] = i;
            L[i]   = low;
          }
          else
          {
            L[i] = qlen;
            qlen++;
            mwm_heap_down(i,n,Q,d,L);
          }//if-

          //updates augment tree
          out[iperm[i]] = k;
          pr [iperm[i]] = j;

        }//for-kk all heap objects

        #ifdef MATCH_DEBUG
        printf("After heapify.  qlen: %d \n",qlen);
        printf("After heapify.  up: %d low: %d \n",
            up, low);
        #endif

        for(Int jdum = 0; jdum < num; jdum++)
        { 
          //If top Q is empty, need to fill
          if(low == up)
          {
            if(qlen == 0)
            {
              //no more nodes
              //goto L1000;
              break;
            }
            i = Q[0];
            if(d[i] >= csp)
            {
              //no more nodex
              //goto L1000;
              break;
            }
            dmin = d[i];

            bool has_more_nodes = true;
            while(has_more_nodes == true)
            {
              mwm_heap_del_root(qlen,n,Q,d,L);
              low--;
              Q[low] = i;
              L[i]   = low;
              if(qlen == 0)
              {
                has_more_nodes = false;
                break;
              }
              i = Q[0];
              if(d[i] > dmin)
              {
                has_more_nodes = false;
              }
            }//while-have more nodes
          }//if(low==up) Q2 empty

          q0 = Q[up-1];
          Entry dq0 = d[q0];

          //if path to q0 is longer than shortest augment path
          //maybe a break ??
          if(dq0 >= csp)
          {
            //goto L1000;
            break;
          }
          up--;

          //scal column that matech with row q0
          j = iperm[q0];
          Entry vj = dq0 - val[jperm[j]] + U[q0];
          for( k = col_ptr[j]; k < col_ptr[j+1]; k++)
          {
            i = row_idx[k];
            //check if located in the Q
            if(L[i] >= up)
            {
              continue;
            }
            //update cost
            Entry dnew_k = vj+ val[k]-U[i];
            //if newcost is more continue
            if(dnew_k >= csp)
            {
              continue;
            }

            //if row is not already matched
            if(iperm[i] == -1)
            {
              csp = dnew_k;
              isp = k;
              jsp = j;
            }
            else
            {
              //check if dnew is smaller
              Entry di = d[i];
              //if not smaller continue
              if(di <= dnew_k)
              {
                continue;
              }
              //if not already in heap (lower heap)
              //it must already be in upper heap (ddmin)
              if(L[i] >= low)
              {
                continue;
              }
              d[i] = dnew_k;
              //if new value is less than min, 
              //needs to be moved to upper heap
              if(dnew_k <= dmin)
              {
                Int lpos = L[i];
                //check that it really is in the heap
                //this check may be overkill
                if(lpos != -1)
                {
                  //extract and reheap at lpos
                  mwm_heap_extract_reheap(lpos,qlen,n,Q,d,L);
                }//if - in Q1
                low--;
                Q[low] = i;
                L[i]   = low;
              }//if = dnew <=dmin, move
              else
              {
                //not in Q already
                if(L[i] != -1)
                {
                  L[i] = qlen;
                  qlen++;
                }
                mwm_heap_down(i,n,Q,d,L);
              }
              out[iperm[i]] = k;
              pr [iperm[i]]  = j;
            }//if(row is or isnot already matched)
          }//for--upated all connecting nodes
        }
      }//if-column not matched


      //find whol augmenting path by going backwards
      //update iperm and jperm
      if(csp != INF)
      {
        num++;
        i          = row_idx[isp];
        iperm[i]   = jsp;
        jperm[jsp] = isp;
        Int j       = jsp;

        for(Int jdum  = 0; jdum < num; jdum++)
        {
          Int jj = pr[j];
          //if at root
          if(jj == -1)
          { 
            break;
          }

          Int k_j   = out[j];
          i         = row_idx[k_j];
          iperm[i]  = jj;
          jperm[jj] = k_j;
          j         = jj;
        }//end for--jdum pathtrace

        //update U for row in Q(up:n)
        //?? < n || <= n
        for(Int kk = up; kk < n; kk++)
        {
          i = Q[kk];
          U[i] = U[i] + d[i] - csp;
        }//end for -- kk update Q(up:)
      }//if csp!=inf

      //wipe clean upper heap
      //?? <n || <= n
      for(Int kk = low; kk <n; kk++)
      {
        i    = Q[kk];
        d[i] = INF;
        L[i] = -1;
      }
      //wipe clean lower heap
      for( k = 0; k < qlen; k++)
      {
        i    = Q[k];
        d[i] = INF;
        L[i] = -1;

      }
    }//for--outer most loop over column matches

    //set dual column variable in d(1:n)
    for(Int j = 0; j < n; j++)
    {
      k = jperm[j];
      if(k != -1)
      {
        d[j] = val[k] - U[row_idx[k]];
      }
      else
      {
        d[j] = (Entry) 0;
      }
      if(iperm[j] == -1)
      {
        U[j] = (Entry) 0;
      }
    }//end for-j , set dual column variables

    //If we where un successful !
    if(num != n)
    {
      //clear j for workspace
      for(Int j = 0; j < n; j++)
      {
        jperm[j] = -1;
      }//for- j, clear j for work space

      k = 0;
      for(i = 0; i < n; i++)
      {
        if(iperm[i] == -1)
        {
          k++;
          out[k] = i;
        }
        else
        {
          Int j = iperm[i];
          jperm[j] = i;
        }
      }//for=i, find row not matched

      k = 0;
      for(Int j=0; j <n ; j++)
      {
        if(jperm[j] != -1)
        {
          continue;
        }
        k++;
        Int jdum = out[k];
        iperm[jdum] = -j;
      }
    }//if-num!=n, we where unsuccessful

    //Done!
    delete [] Q;
    delete [] out;
    return 0;
  }//end mwm_diag_prod()

}//end namespace mwm_order
#endif
