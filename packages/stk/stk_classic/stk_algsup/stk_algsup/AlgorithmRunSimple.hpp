/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_algsup_AlgorithmRunSimple_hpp
#define stk_algsup_AlgorithmRunSimple_hpp

#include <string>

#include <Shards_ArrayVector.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_algsup/AlgorithmRunner.hpp>

namespace stk_classic {
namespace mesh {
/** \ingroup stk_mesh_module
 *  \brief  Field with defined data type and multi-dimensions (if any)
 *  \todo REFACTOR  This convenience class is not appropriately named.
 */
template< class F, class R >
class GatherField
  : public Field<typename FieldTraits<F>::data_type *, R>
{
#ifndef DOXYGEN_COMPILE
       
  public:
    typedef Field<typename FieldTraits<F>::data_type *, R> PtrField;
           
    typedef F RelatedField;
    typedef R Relation;
               
  private:
               
    ~GatherField();
    GatherField();
    GatherField( const GatherField & );
    GatherField & operator = ( const GatherField & );
                       
#endif /* DOXYGEN_COMPILE */
};
} // namespace mesh

struct AlgField
{
  virtual ~AlgField()
  {}
  
  virtual void set(const stk_classic::mesh::Bucket &bucket) = 0;
};


typedef std::vector<AlgField *> AlgFieldVector;


class Algorithm
{
public:
  Algorithm(stk_classic::mesh::MetaData &meta_data)
    : m_metaData(meta_data)
  {}

  virtual ~Algorithm()
  {}

  void add_field(AlgField *field) {
    m_fieldVector.push_back(field);
  }

  virtual void apply( const AlgorithmWork & ) = 0 ;

//   static void number_cruncher(.....);

  
//   void run(const stk_classic::mesh::Bucket& bucket, int begin, int end ) {
//     run(bucket, begin, end);
//   }
  
public:
  const stk_classic::mesh::MetaData &   m_metaData;

private:
  AlgFieldVector                m_fieldVector;
};


template <class T>
struct AlgFieldPtr : public AlgField
{
  typedef typename stk_classic::mesh::FieldTraits<T>::data_type Scalar;
  
  AlgFieldPtr(Algorithm *algorithm, const char *name)
    : m_field(*algorithm->m_metaData.get_field<T>(name))
  {
    algorithm->add_field(this);
  }
  
  AlgFieldPtr(Algorithm *algorithm, const char *name, stk_classic::mesh::FieldState field_state)
    : m_field(algorithm->m_metaData.get_field<T>(name)->field_of_state(field_state))
  {
    algorithm->add_field(this);
  }
  
  AlgFieldPtr(Algorithm *algorithm, T &field)
    : m_field(field)
  {
    algorithm->add_field(this);
  }

  AlgFieldPtr(Algorithm *algorithm, T &field, stk_classic::mesh::FieldState field_state)
    : m_field(field.field_of_state(field_state))
  {
    algorithm->add_field(this);
  }

  virtual ~AlgFieldPtr()
  {}
  
  virtual void set(const stk_classic::mesh::Bucket &bucket) {
    m_ptr = stk_classic::mesh::field_data<T>(m_field, bucket.begin());
  }
  
  T &                   m_field;
  Scalar *              m_ptr;
};


template< class ArrayType >
struct ToArrayVector
{};


template< typename Scalar , shards::ArrayOrder Order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7>
struct ToArrayVector< shards::Array<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> >
{
  typedef shards::ArrayVector<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> ArrayVectorType;
};


template <class T>
struct DimTagInfo;


template <>
struct DimTagInfo<stk_classic::mesh::Cartesian>
{
  static unsigned value(const stk_classic::mesh::Bucket &bucket) {
    return 3;
  }
};


template <>
struct DimTagInfo<void>
{
  static unsigned value(const stk_classic::mesh::Bucket &bucket) {
    return 0;
  }
};


template <class RelationField>
class AlgFieldGather;


template< typename Scalar , 
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7, class Relation>
struct AlgFieldGather< stk_classic::mesh::GatherField <stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, Relation> > : public AlgField
{
  typedef stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> T;
  typedef typename stk_classic::mesh::GatherField<stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, Relation>::PtrField PtrField;
  typedef typename ToArrayVector<typename shards::ArrayAppend< typename shards::ArrayAppend< shards::Array<Scalar,shards::FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, Relation>::type, stk_classic::mesh::EntityDimension >::type>::ArrayVectorType ScratchArray;
  
  AlgFieldGather(Algorithm *algorithm, const char *name)
    : m_field(*algorithm->m_metaData.get_field<PtrField>(name)),
      m_array()
  {
    algorithm->add_field(this);
  }
  
  AlgFieldGather(Algorithm *algorithm, const char *name, stk_classic::mesh::FieldState field_state)
    : m_field(algorithm->m_metaData.get_field<PtrField>(name)->field_of_state(field_state)),
      m_array()
  {
    algorithm->add_field(this);
  }
  
//   AlgFieldGather(Algorithm *algorithm, T &field)
//     : m_field(field),
//       m_array()
//   {
//     algorithm->add_field(this);
//   }

//   AlgFieldGather(Algorithm *algorithm, T &field, stk_classic::mesh::FieldState field_state)
//     : m_field(field.field_of_state(field_state)),
//       m_array()
//   {
//     algorithm->add_field(this);
//   }

  virtual ~AlgFieldGather()
  {}
  
  virtual void set(const stk_classic::mesh::Bucket &bucket) {
    m_begin = stk_classic::mesh::field_data<PtrField>(m_field, bucket.begin());
    m_end = stk_classic::mesh::field_data<PtrField>(m_field, bucket.end());
    
    unsigned dims[8];
    dims[0] = DimTagInfo<Tag1>::value(bucket);
    dims[1] = DimTagInfo<Tag2>::value(bucket);
    dims[2] = DimTagInfo<Tag3>::value(bucket);
    dims[3] = DimTagInfo<Tag4>::value(bucket);
    dims[4] = DimTagInfo<Tag5>::value(bucket);
    dims[5] = DimTagInfo<Tag6>::value(bucket);
    dims[6] = DimTagInfo<Tag7>::value(bucket);
    dims[7] = 0;

    stk_classic::mesh::BucketArray<PtrField> gather_array(m_field, bucket);

    dims[stk_classic::mesh::FieldTraits<T>::Rank] = gather_array.dimension(0);
    dims[stk_classic::mesh::FieldTraits<T>::Rank + 1] = gather_array.dimension(1);

    m_array.resize(&dims[0]);
    m_size = 1;
    for (int i = 0; i < ScratchArray::Rank - 2; ++i)
      m_size *= m_array.dimension(i);
    
    m_ptr = m_array.contiguous_data();
  }

  void fill(const Scalar &value) {
    Scalar *d = m_ptr;
    for (Scalar **p = m_begin; p != m_end; ++p)
      for (Scalar *q = *p; q != *p + m_size; ++q)
        *d++ = value;
  }
  
  void gather() {
    Scalar *d = m_ptr;
    for (Scalar **p = m_begin; p != m_end; ++p)
      for (Scalar *q = *p; q != *p + m_size; ++q)
        *d++ = *q;
  }
  
  void scatter() 
  {
    Scalar *d = m_ptr;
    for (Scalar **p = m_begin; p != m_end; ++p)
      for (Scalar *q = *p; q != *p + m_size; ++q)
        *q = *d++;    
  }

  void assemble() 
  {
    Scalar *d = m_ptr;
    for (Scalar **p = m_begin; p != m_end; ++p)
      for (Scalar *q = *p; q != *p + m_size; ++q)
        *q += *d++;    
  }

  PtrField &            m_field;
  Scalar **             m_begin;
  Scalar **             m_end;
  Scalar *              m_ptr;
  unsigned              m_size;
  ScratchArray          m_array;
};


template <class RelationField>
class AlgFieldX;


template< typename Scalar , 
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7, class Relation>
struct AlgFieldX< stk_classic::mesh::GatherField <stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, Relation> > : public AlgField
{
  typedef stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> T;
  typedef typename stk_classic::mesh::GatherField<stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, Relation>::PtrField PtrField;
  typedef typename shards::ArrayAppend< typename shards::ArrayAppend< shards::Array<Scalar,shards::FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>, Relation>::type, stk_classic::mesh::EntityDimension >::type Array;
  
  AlgFieldX(Algorithm *algorithm, const char *name)
    : m_field(*algorithm->m_metaData.get_field<PtrField>(name)),
      m_begin(0),
      m_end(0),
      m_ptr(0)
  {
    algorithm->add_field(this);
  }
  
  AlgFieldX(Algorithm *algorithm, const char *name, stk_classic::mesh::FieldState field_state)
    : m_field(algorithm->m_metaData.get_field<PtrField>(name)->field_of_state(field_state)),
      m_begin(0),
      m_end(0),
      m_ptr(0)
  {
    algorithm->add_field(this);
  }
  
//   AlgFieldX(Algorithm *algorithm, T &field)
//     : m_field(field),
//       m_array()
//   {
//     algorithm->add_field(this);
//   }

//   AlgFieldX(Algorithm *algorithm, T &field, stk_classic::mesh::FieldState field_state)
//     : m_field(field.field_of_state(field_state)),
//       m_array()
//   {
//     algorithm->add_field(this);
//   }

  virtual ~AlgFieldX()
  {}
  
  virtual void set(const stk_classic::mesh::Bucket &bucket) {
    m_begin = stk_classic::mesh::field_data<PtrField>(m_field, bucket.begin());
    m_end = stk_classic::mesh::field_data<PtrField>(m_field, bucket.end());
    m_ptr = m_begin;
  }

//   PtrField::Array get(unsigned i) 
//   {
//     return PtrField::Array(m_begi
//   void fill(const T::Array &value) {
//     std::fill(m_ptr, m_ptr + m_array.size(), value);
//   }
  
  PtrField &            m_field;
  Scalar **             m_begin;
  Scalar **             m_end;
  Scalar **             m_ptr;
};


template <class T>
struct AlgFieldArray : public AlgField, public stk_classic::mesh::BucketArray<T>
{
  typedef stk_classic::mesh::BucketArray<T> Array;
  
  AlgFieldArray(Algorithm *algorithm, const char *name)
    : m_field(*algorithm->m_metaData.get_field<T>(name))
  {
    algorithm->add_field(this);
  }
  
  AlgFieldArray(Algorithm *algorithm, const char *name, stk_classic::mesh::FieldState field_state)
    : m_field(algorithm->m_metaData.get_field<T>(name)->field_of_state(field_state))
  {
    algorithm->add_field(this);
  }
  
  AlgFieldArray(Algorithm *algorithm, T &field)
    : m_field(field)
  {
    algorithm->add_field(this);
  }

  AlgFieldArray(Algorithm *algorithm, T &field, stk_classic::mesh::FieldState field_state)
    : m_field(field.field_of_state(field_state))
  {
    algorithm->add_field(this);
  }

  virtual void set(const stk_classic::mesh::Bucket &bucket) {
    Array::setup(m_field, bucket);
  }
  
  T &           m_field;
};


template <class Field>
class FillFieldAlgorithm;

template< typename Scalar , 
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7>
class FillFieldAlgorithm< stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
{
public:
  typedef stk_classic::mesh::Field<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> Field;
  typedef shards::ArrayVector<Scalar,shards::FortranOrder, Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> FillArray;
  
  FillFieldAlgorithm(Field &field, const FillArray &fill_value)
    : m_field(field),
      m_fillValue(),
      m_value(0)
  {
    std::vector<unsigned> dims;
    fill_value.dimensions(dims);
    
    m_fillValue.resize(&dims[0]);
    std::copy(fill_value.contiguous_data(), fill_value.contiguous_data() + fill_value.size(), m_fillValue.contiguous_data());
  }

  FillFieldAlgorithm(Field &field, const Scalar &value)
    : m_field(field),
      m_fillValue(),
      m_value(value)
  {}

  enum { chunk_size = 0 }; ///< Don't slice the buckets

  void apply( const stk_classic::AlgorithmWork & work )
  {
    //const stk_classic::mesh::Bucket& bucket = work.bucket ;
    if (m_fillValue.size()) { 
      Scalar *begin_p = stk_classic::mesh::field_data<Field>(m_field, work.bucket_slice_begin);
      Scalar *end_p = stk_classic::mesh::field_data<Field>(m_field, work.bucket_slice_end);
      for (Scalar *p = begin_p; p != end_p; p += m_fillValue.size())
        std::copy(m_fillValue.contiguous_data(), m_fillValue.contiguous_data() + m_fillValue.size(), p);
    }
    else {        
      Scalar *begin_p = stk_classic::mesh::field_data<Field>(m_field, work.bucket_slice_begin);
      Scalar *end_p = stk_classic::mesh::field_data<Field>(m_field, work.bucket_slice_end);
      std::fill(begin_p, end_p, m_value);
    }
  }

private:
  Field &               m_field;
  FillArray             m_fillValue;
  Scalar                m_value;
};

} // namespace stk_classic

#endif
