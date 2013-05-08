/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_TransferP2P_hpp
#define  STK_TransferP2P_hpp

#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/Comm.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

//*********************************************************************
//
//  ........................... NOT FOR PUBLIC USE .....................
//
//  NOTE:
//
//    These routines are used internally by the Transfer class and are 
//    not ment to be used outside of the class. Thus there is no 
//    documentation here on how to call these functions, the
//    documentation is in the source.
//
//    This was done just to make unit testing easier. Otherwise the
//    functions in this namespace STK_TransferP2P would be in an
//    anonomous namespace in the file in which they are used.
//
//    If you want to transfer Point-To-Point call the function
//    Transfer::PointToPoint() in the Transfer class which will
//    then call these functions in the correct order and wtih the
//    coorrect data.
//
//*********************************************************************
//
//
// Local Typedefs.  These might need to be lifted to be based on template
// parameters.
typedef stk::search::ident::IdentProc<uint64_t,unsigned>            IdentProc;
typedef std::vector<std::pair<IdentProc, IdentProc> >               IdentProcRelation;
typedef stk::search::box::SphereBoundingBox<IdentProc,float,3>      BoundingBox;
typedef std::map<unsigned, std::vector<double> >                    PointMap;

namespace STK_TransferP2P {

typedef Intrepid::FieldContainer<double> MDArray;


int LU_decomp(double A[9], int piv[3], int* sign);
int LU_solve (const double A[9], const int piv[3], double b[3]);
std::vector<double> solve_3_by_3_with_LU(const MDArray             M, 
                                         const std::vector<double> x);

template <unsigned DIM>
void point_to_point_coarse_search(IdentProcRelation &RangeToDomain,
                                  const PointMap    &ToPoints,
                                  const MDArray     &FromPoints,
                                  const BoundingBox::Data radius,
                                  const stk::ParallelMachine comm);
template <unsigned DIM>
void linear_interpolation (MDArray &ToValues,
                    const MDArray &FromValues,
                    const IdentProcRelation &RangeToDomain,
                    const MDArray &ToPoints,
                    const MDArray &FromPoints,
                    const stk::ParallelMachine  comm);
  
template <unsigned DIM>
void filter_with_fine_search(IdentProcRelation &RangeToDomain,
                                const PointMap &ToPoints,
                                const MDArray &FromPoints,
                                const stk::ParallelMachine  comm) ;

template <unsigned DIM>
void delete_range_points_found(PointMap &ToPoints,
                               const IdentProcRelation &RangeToDomain,
                               const stk::ParallelMachine  comm) ;

template <unsigned DIM>
void convert_to_map(PointMap      &map_points,
                    const MDArray &Points) ;

}
#endif

