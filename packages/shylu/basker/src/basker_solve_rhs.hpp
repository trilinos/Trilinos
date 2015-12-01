#ifndef BASKER_SOLVE_RHS_HPP
#define BASKER_SOLVE_RHS_HPP

/*Basker Includes*/
#include "basker_decl.hpp"
#include "basker_matrix_decl.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_types.hpp"
#include "basker_util.hpp"

/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif

/*System Includes*/
#include <iostream>
#include <string>

//#define BASKER_DEBUG_SOLVE_RHS

using namespace std;

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::test_solve()
  {
    ENTRY_1DARRAY  x_known;
    ENTRY_1DARRAY  x;
    ENTRY_1DARRAY  y;


    printf("test_solve called \n");

    BASKER_ASSERT(gn > 0, "solve testsolve gn");
    MALLOC_ENTRY_1DARRAY(x_known, gn);
    init_value(x_known, gn , (Entry)1.0);


    //temp
    for(Int i = 0; i < gn; i++)
      {
	//x_known(i) = (Entry)(i+1);
        x_known(i) = (Entry) 1.0;
      }
    //JDB: used for other test
    //permute(x_known, order_csym_array, gn);



    MALLOC_ENTRY_1DARRAY(x, gn);
    init_value(x, gn, (Entry) 0.0);
    BASKER_ASSERT(gm > 0, "solve testsolve gm");
    MALLOC_ENTRY_1DARRAY(y, gm);
    init_value(y, gm, (Entry) 0.0);
    
    if(btf_nblks > 0)
      {
	sort_matrix(BTF_C);
	//printMTX("C_BEFORE_SOLVE.mtx", BTF_C);
      }

    if(Options.btf == BASKER_TRUE)
      {
      
	printf("btf_tabs_offset: %d ", btf_tabs_offset);
        printf("btf_nblks: %d \n", btf_nblks);
	if(btf_tabs_offset != 0)
	  {
            printf("BTF_A spmv\n");
	    spmv(BTF_A, x_known,y);
            if(btf_nblks> 1)
              {
                printf("btf_B spmv \n");
                spmv(BTF_B, x_known, y);
              }
          }
        if(btf_nblks > 1)
          {
            
            printf("btf_c spmv \n");
            spmv(BTF_C, x_known, y);
          }
	//return -1;
      }
    else
      {
        //printf("other\n");
	//spmv(BTF_A, x_known,y);
      }
    
    //return -1;

    printf("DEBUG\n");
    
    printf("\n Before Test Points \n");
    printf("i: %d x: %f y: %f \n", 0, x_known(0), y(0));
    printf("i: %d x: %f y: %f \n", 24, x_known(24), y(24));
    
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("Known Solution: \n");
    for(Int i = 0; i < gn; i++)
      {
	printf("%f, " , x_known(i));
      }
    printf("\n\n");
    printf("RHS: \n");
    for(Int i =0; i < gm; i++)
      {
	printf("%f, ", y(i)); 
      }
    printf("\n\n");
    #endif


    if(Options.btf == BASKER_FALSE)
      {
	if(btf_tabs_offset != 0)
	  {
	    serial_solve(y,x);
	  }
	//printf("After serial solve\n");
	//printf("i: %d x: %f y: %f \n", 0, x(0), y(0));
	//printf("i: %d x: %f y: %f \n", 24, x(24), y(24));
      }
    else
      {
	//A\y -> y
	//serial_btf_solve(y,x);
	serial_btf_solve(y,x);

	//printf("After btf solve\n");
	//  printf("i: %d x: %f y: %f \n", 0, x(0), y(0));
	//  printf("i: %d x: %f y: %f \n", 24, x(24), y(24));
      }


    Entry diff =0.0;
  
    for(Int i = 0; i < gn; i++)
      {
	diff += (x_known(i) - x(i));
      }
    diff = diff/(Entry) gn;

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("Solve Compare: \n");
    for(Int i = 0; i < gn; i++)
      {
	printf("%d %f %f \n", 
	       i, x_known(i), x(i));
 
      }
    printf("\n\n");
    #endif



    printf("\n Test Points \n");
    printf("i: %d x: %f %f \n", 0, x_known(0), x(0));
    printf("i: %d x: %f %f \n", 10, x_known(10), x(10));
    printf("i: %d x: %f %f \n", 24, x_known(24), x(24));
    printf("\n");
    printf("TEST_SOLVE: ||x-x||/||x| = %e", diff);
    printf("\n");

    if((diff > -1e-2) && (diff < 1e-2))
      {
        printf("TEST PASSED \n");
      }  

    return 0;
  }//end test_solve


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interface
  (
   Entry *_x, //Solution (len = gn)
   Entry *_y
   )
  {
  
    //===== Move to view
    ENTRY_1DARRAY  x;
    ENTRY_1DARRAY  y;

    MALLOC_ENTRY_1DARRAY(x, gn);
    MALLOC_ENTRY_1DARRAY(y, gm);
    
    for(Int i =0; i < gn; i++)
      {
	x(i) = (Entry) 0;
	y(i) = (Entry) _y[i];
      }

    //===== Permute
    //printf("Permute RHS\n");
    //==== Need to make this into one global perm
    if(match_flag == BASKER_TRUE)
      {
	//printf("match order\n");
	//printVec(order_match_array, gn);
	permute_inv(y,order_match_array, gn);
      }
    if(btf_flag == BASKER_TRUE)
      {
	//printf("btf order\n");
	//printVec(order_btf_array, gn);
	permute_inv(y,order_btf_array, gn);
      }
    if(nd_flag == BASKER_TRUE)
      {
	//printf("ND order \n");
	//printVec(part_tree.permtab, gn);
	permute_inv(y,part_tree.permtab, gn);
      }
    if(amd_flag == BASKER_TRUE)
      {
	//printf("AMD order \n");
	//printVec(order_csym_array, gn);
	//FILE *fpamd;
	//fpamd = fopen("amd.csv", "w");
	//for(Int i = 0; i < gn; i++)
	// {
	//    fprintf(fpamd, "%d \n", 
	//	    order_csym_array(i));
	//  }
	//fclose(fpamd);
	    
	permute_inv(y,order_csym_array, gn);
      }

    solve_interface(x,y);

    //Inverse perm
    if(match_flag == BASKER_TRUE)
      {
	//printf("match order\n");
	//printVec(order_match_array, gn);
	permute(x,order_match_array, gn);
      }
    if(btf_flag == BASKER_TRUE)
      {
	//printf("btf order\n");
	//printVec(order_btf_array, gn);
	permute(x,order_btf_array, gn);
      }
    if(nd_flag == BASKER_TRUE)
      {
	//printf("ND order \n");
	//printVec(part_tree.permtab, gn);
	permute(x,part_tree.permtab, gn);
      }
    if(amd_flag == BASKER_TRUE)
      {
	//printf("AMD order \n");
	//printVec(order_csym_array, gn);
	//FILE *fpamd;
	//fpamd = fopen("amd.csv", "w");
	//for(Int i = 0; i < gn; i++)
	// {
	//   fprintf(fpamd, "%d \n", 
	//	    order_csym_array(i));
	//  }
	//fclose(fpamd);
	    
	permute(x,order_csym_array, gn);
      }
   


    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("X: \n");
    for(Int i = 0; i < gn; i++)
      {
	printf("%f, " , x(i));
      }
    printf("\n\n");
    printf("RHS: \n");
    for(Int i =0; i < gm; i++)
      {
	printf("%f, ", y(i)); 
      }
    printf("\n\n");
    #endif


    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interface
  (
   ENTRY_1DARRAY x,//Solution (len = gn)
   ENTRY_1DARRAY y //RHS      (len = gn)
   )
  {
   
    
    //printf("\n Before Test Points \n");
    //printf("i: %d x: %f y: %f \n", 0, x(0), y(0));
    //printf("i: %d x: %f y: %f \n", 24, x(24), y(24));
    
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("X: \n");
    for(Int i = 0; i < gn; i++)
      {
	printf("%f, " , x(i));
      }
    printf("\n\n");
    printf("RHS: \n");
    for(Int i =0; i < gm; i++)
      {
	printf("%f, ", y(i)); 
      }
    printf("\n\n");
    #endif


    if(Options.btf == BASKER_FALSE)
      {
	if(btf_tabs_offset != 0)
	  {
	   
	    serial_solve(y,x);
	    
	    //printf("After serial solve\n");
	    //printf("i: %d x: %f y: %f \n", 0, x(0), y(0));
	    //printf("i: %d x: %f y: %f \n", 24, x(24), y(24));
   
	  }
      }
    else
      {
	//A\y -> y
	//serial_btf_solve(y,x);
	serial_btf_solve(y,x);

	//printf("After btf solve\n");
	// printf("i: %d x: %f y: %f \n", 0, x(0), y(0));
	// printf("i: %d x: %f y: %f \n", 24, x(24), y(24));
   
      }


    //printf("\n After Test Points \n");
    //printf("i: %d x: %f y: %f \n", 0, x(0), y(0));
    //printf("i: %d x: %f y: %f \n", 24, x(24), y(24));
    
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("X: \n");
    for(Int i = 0; i < gn; i++)
      {
	printf("%f, " , x(i));
      }
    printf("\n\n");
    printf("RHS: \n");
    for(Int i =0; i < gm; i++)
      {
	printf("%f, ", y(i)); 
      }
    printf("\n\n");
    #endif

    return 0;
  }
  

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_solve
  (
   ENTRY_1DARRAY y,
   ENTRY_1DARRAY x
   )
  {
    for(Int i =0; i < gn; ++i)
      {
	//x[i] = y[i];
	x(i) = y(i);
	//y[i] = 0.0;
	y(i) = (Entry) 0.0;
      }
    //L\x -> y
    serial_forward_solve(x,y);
    for(Int i =0; i<gn; ++i)
      {
	//x[i] = 0;
	x(i) = 0;
      }
    //U\y -> x
    serial_backward_solve(y,x);

    return 0;
  }//end serial solve()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_btf_solve
  (
   ENTRY_1DARRAY y,
   ENTRY_1DARRAY x
   )
  {

    
    for(Int i = 0; i < gn; ++i)
      {
	x(i) = y(i);
	y(i) = (Entry) 0.0;
      }
    

    //Start in C and go backwards
    //In first level, only due U\L\x->y
    for(Int b = (btf_nblks-btf_tabs_offset)-1;
	b>= 0; b--)
      {

	//---Lower solve
	BASKER_MATRIX &LC = LBTF(b);
	//L\x -> y 
	lower_tri_solve(LC,x,y);

	BASKER_MATRIX &UC = UBTF(b);
	//U\x -> y
	upper_tri_solve(UC,x,y);

       
	//-----Update
	//if(b > btf_tabs_offset)
	  {
	//x = BTF_C*y;
            //printf("spmv tab: %d \n", b+btf_tabs_offset);
	spmv_BTF(b+btf_tabs_offset,
		 BTF_C, y, x);
	  }

	//BASKER_MATRIX &UC = UBTF[b];
	//U\x -> y
	//upper_tri_solve(UC,x,y);

      }



    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, BTF-C Solve \n");
    printf("\n x \n");
    printVec(x, gn);
    printf("\n y \n");
    printVec(y, gn);
    printf("\n\n");
    #endif

    //Update B
    //BTF_B*y -> x
    if(btf_tabs_offset !=  0)
      {
    neg_spmv(BTF_B,y,x);
      }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, SPMV BTF_B UPDATE \n");
    printf("\n x \n");
    printVec(x, gn);
    printf("\n y \n");
    printVec(y, gn);
    printf("\n\n");
    #endif

    //now do the forward backwared solve
    //L\x ->y
    serial_forward_solve(x,y);

    //U\y->x
    serial_backward_solve(y,x);

    //copy lower part down
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("copying lower starting: %d \n",
	   btf_tabs[btf_tabs_offset]);
    #endif
    for(Int i = btf_tabs(btf_tabs_offset); i < gn; ++i)
      {
	//x[i] = y[i];
	x(i) = y(i);
      }
    
    //Comeback and fix
    return 0;
    
  }//end serial_btf_solve


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_forward_solve
  (
   ENTRY_1DARRAY y, 
   ENTRY_1DARRAY x
   )
  {
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Called serial forward solve \n");
    #endif

   
    //Forward solve on A
    for(Int b = 0; b < tree.nblks; ++b)
      {
	#ifdef BASKER_DEBUG_SOLVE_RHS
	printf("Lower Solve blk: %d \n",b);
	#endif

	BASKER_MATRIX &L = LL(b)(0);
	//L\y -> x
	lower_tri_solve(L, y, x);

	//Update offdiag
	for(Int bb = 1; bb < LL_size(b); ++bb)
	  {
	    #ifdef BASKER_DEBUG_SOLVE_RHS
	    printf("Lower Solver Update blk: %d %d \n",
		   b, bb);
	    #endif

	    BASKER_MATRIX &LD = LL(b)(bb);
	    //y = LD*x;
	    neg_spmv(LD, x, y);
	  }
      }
   
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done forward solve A \n");
    printVec(x, gn);
    #endif
    
    return 0;
  }//end serial_forward_solve()

  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_backward_solve
  (
   ENTRY_1DARRAY y,
   ENTRY_1DARRAY x
   )
  {
    
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("called serial backward solve \n");
    #endif


    for(Int b = tree.nblks-1; b >=0; b--)
      {
	//printf("--HERE-\n");
        #ifdef BASKER_DEBUG_SOLVE_RHS
	printf("Upper solve blk: %d \n", b);
	#endif
	
	BASKER_MATRIX &U = LU(b)(LU_size(b)-1);

	//printf("\n--------U------------\n");
	//U.print();

	//U\y -> x
	upper_tri_solve(U,y,x);
	
       
	for(Int bb = LU_size(b)-2; bb >= 0; bb--)
	  {
	    #ifdef BASKER_DEBUG_SOLVE_RHS
	    printf("Upper solver spmv: %d %d \n",
		   b, bb);
	    #endif
	    
	    BASKER_MATRIX &UB = LU(b)(bb);
	    //y = UB*x;
	    neg_spmv(UB,x,y);
	  }
	
      }//end over all blks
    
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done with Upper Solve: \n");
    printVec(x, gn);
    #endif
    
    return 0;
  }//end serial_backward_solve()


  //Horrible, cheap spmv
  //y = M*x
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::spmv(
					BASKER_MATRIX &M,
					ENTRY_1DARRAY x,
					ENTRY_1DARRAY y)
  {
    //Add checks
    //#ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. scol: %d ncol: %d nnz: %d \n",
	   M.scol, M.ncol, M.nnz);
    M.info();
    //#endif

    const Int bcol = M.scol;
    const Int brow = M.srow;
    //for(Int k=M.scol; k < (M.scol+M.ncol); k++)
    for(Int k = 0; k < M.ncol; ++k)
      {
	//printf("k: %d \n", k);
	for(Int i = M.col_ptr(k); i<M.col_ptr(k+1); ++i)
	  {
	    const Int j = M.row_idx(i);
	    //printf("j: %d i: %d idx1: %d idx2: %d \n",
	    //	   j, i, j+brow, k+bcol);
	    
	    //y[j] += M.val[i]*x[k];
	    y(j+brow) += M.val(i)*x(k+bcol);

	  }
      }
    return 0;
  }//spmv

  //Horrible, cheap spmv
  //y = M*x
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::neg_spmv(
					BASKER_MATRIX &M,
					ENTRY_1DARRAY x,
					ENTRY_1DARRAY y)
  {
    //Add checks
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. scol: %d ncol: %d \n",
	   M.scol, M.ncol);
    #endif

    const Int bcol = M.scol;
    const Int brow = M.srow;
    //for(Int k=M.scol; k < (M.scol+M.ncol); k++)
    for(Int k=0; k < M.ncol; ++k)
      {
	//for(Int i = M.col_ptr[k-bcol];
	//  i < M.col_ptr[k-bcol+1]; i++)
	for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
	  {
	    //Int j = M.row_idx[i];
	    const Int j = M.row_idx(i);
	    
	    //y[j] -= M.val[i]*x[k];
	    y(j+brow) -= M.val(i)*x(k+bcol);

	  }
      }

  }//neg_spmv


  //M\x = y
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::lower_tri_solve
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY x, 
   ENTRY_1DARRAY y
   )
  {
    const Int bcol = M.scol;
    const Int brow = M.scol;
    
    //printf("Lower-Tri-Solve-Test, [%d %d %d %d] \n",
    //	  M.srow, M.nrow, M.scol, M.ncol);

    for(Int k = 0; k < M.ncol; ++k)
      {
	//Test if zero pivot value
	#ifdef BASKER_DEBUG_SOLVE_RHS
	BASKER_ASSERT(M.val[M.col_ptr[k]]!=0.0, "LOWER PIVOT 0");
	#endif
	
	if(M.val[M.col_ptr[k]] == 0.0)
	  {
	    printf("Lower Pivot: %d %f \n", 
		   M.row_idx[M.col_ptr[k]],
		   M.val[M.col_ptr[k]]);
	    return -1;
	  }

	//printf("Lower tri.  k: %d out: %f in: %f piv: %f \n",
	//   k+bcol, y[k+bcol], x[k+bcol], M.val[M.col_ptr[k]]);

	//Replace with Entry divide in future
	//y[k+bcol] = x[k+bcol] / M.val[M.col_ptr[k]];
	y(k+brow) = x(k+bcol) / M.val(M.col_ptr(k));
	
	//for(Int i = M.col_ptr[k]+1; i < M.col_ptr[k+1]; i++)
	for(Int i = M.col_ptr(k)+1; i < M.col_ptr(k+1); ++i)
	  {
	   

	    //Int j = gperm[M.row_idx[i]];
	    const Int j = gperm(M.row_idx(i)+brow);
	    
	    #ifdef BASKER_DEBUG_SOLVE_RHS
	    BASKER_ASSERT(j != BASKER_MAX_IDX,"Using nonperm\n");
	    #endif
	    
	    //x[j] -= M.val[i]*y[k+bcol];
	    x(j) -= M.val(i)*y(k+bcol);

	  }//over all nnz in a column

      }//over each column

    return 0;

  }//end lower_tri_solve


  //U\x = y
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::upper_tri_solve
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY x,
   ENTRY_1DARRAY y
   )
  {
    const Int bcol = M.scol;
    const Int brow = M.srow;
    
    //printf("Upper Tri Solve, scol: %d ncol: %d \n",
    //	   M.scol, M.ncol);

    //end over all columns
    for(Int k = M.ncol; k >= 1; k--)
      {
	
	//printf("Upper Tri Solve: k: %d \n", k);

	#ifdef BASKER_DEBUG_SOLVE_RHS
	BASKER_ASSERT(M.val[M.col_ptr[k]-1]!=0.0,"UpperPivot\n");
	#endif

	//if(M.val[M.col_ptr[k]-1]==0.0)
	if(M.val(M.col_ptr(k)-1)==0)
	  {
	    printf("Upper pivot: %d %f \n",
		   M.row_idx[M.col_ptr[k]-1],
		   M.val[M.col_ptr[k]-1]);
	    return -1;
	  }


	//printf("TEST, k: %d  out: %f in: %f pivot: %f\n", 
	//   k, y[k+bcol-1], x[k+bcol-1], 
	//   M.val[M.col_ptr[k]-1]);

	//Comeback and do with and entry divide
	//y[k+bcol-1] = x[k+bcol-1] / M.val[M.col_ptr[k]-1];
	y(k+brow-1)  =  x(k+bcol-1) / M.val(M.col_ptr(k)-1);
	
	
	//for(Int i = M.col_ptr[k]-2; i >= M.col_ptr[k-1]; i--)
	for(Int i = M.col_ptr(k)-2; i >= M.col_ptr(k-1); --i)
	  {
	    //Int j = M.row_idx[i];
	    const Int j = M.row_idx(i);

	    // printf("Updating row_idx: %d %f %f \n",
	    //	   j, x[j], M.val[i]*y[k+bcol-1]);

	    //x[j] -=  M.val[i]*y[k+bcol-1];
	    x(j+brow) -= M.val(i) * y(k+bcol-1);
	  }

      }//end over all columns

    return 0;
  }//end upper_tri_solve

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::spmv_BTF
  (
   Int tab,
   BASKER_MATRIX &M,
   ENTRY_1DARRAY x,
   ENTRY_1DARRAY y
   )
  {
    //Tab = block in    
    const Int bcol = btf_tabs(tab)- M.scol;
    const Int brow = M.srow;
    const Int ecol = btf_tabs(tab+1) - M.scol;
    const Int erow = btf_tabs(tab-1);

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("BTF_UPDATE, TAB: %d [%d %d] [%d %d] \n",
	   tab, brow, erow, bcol, ecol);
    #endif

    //loop over each column
    for(Int k = bcol; k < ecol; ++k)
      {
	//for(Int i = M.col_ptr[k]; i < M.col_ptr[k+1]; i++)
	for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
	  {
	    //Int j = M.row_idx[i];
	    const Int j = M.row_idx(i);
	    if(j > erow)
	      {
		#ifdef BASKER_DEBUG_SOLVE_RHS
		printf("break, k: %d j: %d erow: %d\n",
		       k, j, erow);
		#endif
		//break; //breaks for 1 colummn
		continue;
	      }

	    #ifdef BASKER_DEBUG_SOLVE_RHS
	    printf("BTF_UPDATE-val, j: %d y: %f x: %f, val: %f \n",
		   j, y[j], x[k+M.scol], M.val[i]);
	    #endif
	    //for now just do a single function with zero
	    //y[j] -= M.val[i]*x[k+M.scol];

	    y(j+brow) -= M.val(i)*x(k+M.scol);
	  }//over all nnz in row
      }
    return 0;
  }//end spmv_BTF();

  
 
}//end namespace BaskerNS
#endif //end ifndef basker_solver_rhs
