/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef USECASE_14_HPP
#define USECASE_14_HPP

int use_case_14();

#include <app/UseCase_14_Common.hpp>
#include <app/UseCase_14_Fields.hpp>

typedef stk::mesh::Field<double, stk::mesh::Cartesian>       VectorField ;

class MyHexInternalForceAlg
{
public:
  MyHexInternalForceAlg(const stk::app::Fields &fields,
			lame::matParams &matParameters,
			lame::Material *matModel,
                        stk::mesh::MetaData &md);
 
  enum { maximum_entity_count = 1000 };

  void apply( stk::mesh::Bucket::iterator ,
              stk::mesh::Bucket::iterator ,
              const stk::mesh::PartVector & ) const ;

private:
  MyHexInternalForceAlg& operator=(const MyHexInternalForceAlg&);
  MyHexInternalForceAlg(const MyHexInternalForceAlg&);

  APSHex8ug                             hex_element;
  lame::matParams &                     materialParameters;
  stk::mesh::Property<double> *         delta_t;
  stk::mesh::Property<lame::MatProps> * m_materialProperties;
  lame::Material *                      matmodel;

  stk::app::Fields                      m_fields;
};

//--------------------------------------------------------------------
//--------------------------------------------------------------------

class MyNodalForceScatterAlg
{
public:
  MyNodalForceScatterAlg(const stk::app::Fields &fields, stk::mesh::EntityRank element_rank);

  enum { maximum_entity_count = 0 }; /**< Don't slice the buckets */

  void apply( stk::mesh::Bucket::iterator ,
              stk::mesh::Bucket::iterator ) const ;


private:
  MyNodalForceScatterAlg& operator=(const MyNodalForceScatterAlg&);
  MyNodalForceScatterAlg(const MyNodalForceScatterAlg&);

  ElementNodeVectorField *      force_new_field;
  VectorField *                 fint_field;
  stk::mesh::EntityRank         m_elementRank;
};

#undef INLINE /* */

#endif

