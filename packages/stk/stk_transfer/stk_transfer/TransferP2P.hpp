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
#include <stk_mesh/base/EntityKey.hpp>


#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

#include <boost/smart_ptr/shared_ptr.hpp>

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
typedef stk::search::ident::IdentProc<stk::mesh::EntityKey,unsigned> IdentProc;
typedef std::vector<std::pair<IdentProc, IdentProc> >                IdentProcRelation;
typedef stk::search::box::SphereBoundingBox<IdentProc,float,3>       BoundingBox;
typedef std::map<unsigned, std::vector<double> >                     PointMap;

namespace STK_TransferP2P {


typedef boost::shared_ptr<stk::mesh::BulkData> BulkDataPtr;
typedef boost::shared_ptr<stk::mesh::MetaData> MetaDataPtr;
typedef stk::mesh::Field<double,stk::mesh::Cartesian>    VectorFieldType ;
typedef std::vector<stk::mesh::Entity> EntityVec;


class STKMesh  {
public :
  typedef EntityVec::const_iterator iterator;
  STKMesh(MetaDataPtr meta,
          BulkDataPtr bulk,
          EntityVec   &ent,
          VectorFieldType &coord,
          VectorFieldType &val);
  STKMesh(const STKMesh &M);
  ~STKMesh();
  BulkDataPtr &BulkData();
  iterator begin() const;
  iterator   end() const;
  const double *Coord(const iterator i) const ;
  const double *Coord(const stk::mesh::EntityKey i) const;
  const double *Value(const stk::mesh::EntityKey i) const;
  stk::mesh::EntityKey Key(iterator i) const;
  VectorFieldType &Coord();
  VectorFieldType &Value();
private :
  STKMesh ();
  STKMesh &operator=(const STKMesh&);
  MetaDataPtr meta_data;
  BulkDataPtr bulk_data;
  EntityVec   entities;
  VectorFieldType &coordinates_field;
  VectorFieldType &values_field;
};



typedef Intrepid::FieldContainer<double> MDArray;


int LU_decomp(double A[9], int piv[3], int* sign);
int LU_solve (const double A[9], const int piv[3], double b[3]);
std::vector<double> solve_3_by_3_with_LU(const MDArray             M, 
                                         const std::vector<double> x);

template <unsigned DIM, class PointData>
void point_to_point_coarse_search(IdentProcRelation &RangeToDomain,
                                  const PointMap    &ToPoints,
                                  const PointData   &FromPoints,
                                  const BoundingBox::Data radius,
                                  const stk::ParallelMachine comm);
template <unsigned DIM, class MeshClass>
void linear_interpolation (MDArray          &ToValues,
                    const MeshClass         &FromValues,
                    const IdentProcRelation &RangeToDomain,
                    const MDArray           &ToPoints,
                    const MeshClass         &FromPoints,
                    const stk::ParallelMachine  comm);
  
template <unsigned DIM, class MeshClass>
void filter_with_fine_search(IdentProcRelation &RangeToDomain,
                                const PointMap &ToPoints,
                                const MeshClass &FromPoints,
                                const stk::ParallelMachine  comm) ;

template <unsigned DIM>
void delete_range_points_found(PointMap &ToPoints,
                               const IdentProcRelation &RangeToDomain,
                               const stk::ParallelMachine  comm) ;

template <unsigned DIM>
void convert_to_map(PointMap      &map_points,
                    const MDArray &Points) ;

STKMesh convert_points_to_mesh(const MDArray &Coords, 
                               const MDArray &Values,
                               const stk::ParallelMachine  comm);

void copy_domain_to_range_processors(STKMesh                   &Mesh,
                                     const IdentProcRelation   &RangeToDomain,
                                     const std::string         &transfer_name,
                                     const stk::ParallelMachine comm);
}
#endif

