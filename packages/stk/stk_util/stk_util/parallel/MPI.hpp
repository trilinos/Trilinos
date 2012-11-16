#ifndef STK_UTIL_PARALLEL_MPI_hpp
#define STK_UTIL_PARALLEL_MPI_hpp

#include <stk_util/stk_config.h>
#if defined( STK_HAS_MPI )

#include <mpi.h>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <complex>

namespace sierra {
namespace MPI {

template<class T>
inline
void
  mpi_real_complex_sum(
    void *		invec,
    void *		inoutvec,
    int *		len,
    MPI_Datatype *	datatype)
  {
    std::complex<T> *complex_in = static_cast<std::complex<T> *>(invec);
    std::complex<T> *complex_inout = static_cast<std::complex<T> *>(inoutvec);

    for (int i = 0; i < *len; ++i)
      complex_inout[i] += complex_in[i];
  }

///
/// @addtogroup MPIDetail
/// @{
///


/**
 * @brief Function <code>float_complex_type</code> returns an MPI complex data type for
 * C++.
 *
 * @return	a <code>MPI_Datatype</code> value of the C++ complex MPI data type.
 */
MPI_Datatype float_complex_type();

/**
 * @brief Function <code>double_complex_type</code> returns an MPI complex data type for
 * C++.
 *
 * @return	a <code>MPI_Datatype</code> value of the C++ complex MPI data type.
 */
MPI_Datatype double_complex_type();

/**
 * @brief Function <code>double_complex_sum_op</code> returns a sum operation for the C++
 * complex MPI data type.
 *
 * @return a <code>MPI_Op</code> ...
 */
MPI_Op double_complex_sum_op();

/**
 * @brief Function <code>real_complex_sum_op</code> returns a sum operation for the C++
 * complex MPI data type.
 *
 * @return a <code>MPI_Op</code> ...
 */
template<class T>
inline
MPI_Op real_complex_sum_op()
{
  static MPI_Op s_mpi_real_complex_sum;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Op_create(mpi_real_complex_sum<T>, true, &s_mpi_real_complex_sum);
  }
  return s_mpi_real_complex_sum;
}

/**
 * @brief Member function <b>double_double_int_type</b> ...
 *
 * @return			a <b>MPI_Datatype</b> ...
 */
MPI_Datatype long_long_int_int_type();


/**
 * @brief Member function <b>double_double_int_type</b> ...
 *
 * @return			a <b>MPI_Datatype</b> ...
 */
MPI_Datatype double_double_int_type();


/**
 * @brief Template class <code>loc</code> implements the data structure for the MINLOC and
 * MAXLOC data types.
 *
 */
template <typename T>
struct Loc
{
  Loc()
    : m_value(),
      m_loc(0)
  {}
  
  Loc(const T &value, int loc)
    : m_value(value),
      m_loc(loc)
  {}
  
  T		m_value;
  int		m_loc;
};

struct TempLoc
{
  TempLoc()
    : m_value(),
      m_other(),
      m_loc(0)
  {}
  
  TempLoc(double value, double other, int loc)
    : m_value(value),
      m_other(other),
      m_loc(loc)
  {}
  
  double        m_value;
  double        m_other;
  int		m_loc;
};


/**
 * @brief Traits class <code>Datatype</code> implements a traits class containing two
 * static member functions which return the appropriate MPI data type for the C++ data
 * type.
 *
 * The <b>type()</b> function returns the MPI data type.
 *
 */
template <typename T>
struct Datatype;

template <>
struct Datatype<char>
{
  static MPI_Datatype type() {
    return MPI_CHAR;
  }
};

template <>
struct Datatype<signed char>
{
  static MPI_Datatype type() {
    return MPI_CHAR;
  }
};

template <>
struct Datatype<unsigned char>
{
  static MPI_Datatype type() {
    return MPI_BYTE;
  }
};

template <>
struct Datatype<int>
{
  static MPI_Datatype type() {
    return MPI_INT;
  }
};

template <>
struct Datatype<unsigned int>
{
  static MPI_Datatype type() {
    return MPI_UNSIGNED;
  }
};

template <>
struct Datatype<short>
{
  static MPI_Datatype type() {
    return MPI_SHORT;
  }
};

template <>
struct Datatype<unsigned short>
{
  static MPI_Datatype type() {
    return MPI_UNSIGNED_SHORT;
  }
};

template <>
struct Datatype<long>
{
  static MPI_Datatype type() {
    return MPI_LONG;
  }
};

template <>
struct Datatype<unsigned long>
{
  static MPI_Datatype type() {
    return MPI_UNSIGNED_LONG;
  }
};

// #ifdef MPI_LONG_LONG_INT
// template <>
// struct Datatype<long long>
// {
//   static MPI_Datatype type() {
//     return MPI_LONG_LONG_INT;
//   }
// };
// #endif

// #ifdef MPI_UNSIGNED_LONG_LONG_INT
// template <>
// struct Datatype<unsigned long long>
// {
//   static MPI_Datatype type() {
//     return MPI_UNSIGNED_LONG_LONG_INT;
//   }
// };
// #endif

template <>
struct Datatype<float>
{
  static MPI_Datatype type() {
    return MPI_FLOAT;
  }
};

template <>
struct Datatype<double>
{
  static MPI_Datatype type() {
    return MPI_DOUBLE;
  }
};


template <>
struct Datatype<std::complex<float> >
{
  static MPI_Datatype type() {
    return float_complex_type();
  }
};

template <>
struct Datatype<std::complex<double> >
{
  static MPI_Datatype type() {
    return double_complex_type();
  }
};


template <>
struct Datatype<Loc<int> >
{
  static MPI_Datatype type() {
    return MPI_2INT;
  }
};

template <>
struct Datatype<Loc<short> >
{
  static MPI_Datatype type() {
    return MPI_SHORT_INT;
  }
};

template <>
struct Datatype<Loc<long> >
{
  static MPI_Datatype type() {
    return MPI_LONG_INT;
  }
};

template <>
struct Datatype<Loc<unsigned long> >
{
  static MPI_Datatype type() {
    return MPI_LONG_INT;
  }
};

// #ifdef MPI_LONG_LONG_INT
// template <>
// struct Datatype<Loc<long long> >
// {
//   static MPI_Datatype type() {
//     return long_long_int_int_type();
//   }
// };
// #endif

template <>
struct Datatype<Loc<float> >
{
  static MPI_Datatype type() {
    return MPI_FLOAT_INT;
  }
};

template <>
struct Datatype<Loc<double> >
{
  static MPI_Datatype type() {
    return MPI_DOUBLE_INT;
  }
};

template <>
struct Datatype<TempLoc>
{
  static MPI_Datatype type() {
    return double_double_int_type();
  }
};


/**
 * @brief Function <code>AllReduce</code> copies the source/destination array into a
 * temporary vector and then executed the MPI operation using the temporary as the source.
 *
 * @param mpi_comm	a <code>MPI_Comm</code> value of the MPI communicator.
 *
 * @param op		a <code>MPI_Op</code> value of the MPI operation.
 *
 * @param src_dest	a <code>T</code> pointer to an array to be copied for the source
 *			and used as the destination.
 *
 * @param size		a <code>size_t</code> value of the length of the array pointed to
 *			by <b>src_dest</b>
 *
 */
template<class T>
inline void
AllReduce(MPI_Comm mpi_comm, MPI_Op op, T *src_dest, size_t size)
{
  std::vector<T> source(src_dest, src_dest + size);

  if (MPI_Allreduce(&source[0], &src_dest[0], (int) size, Datatype<T>::type(), op, mpi_comm) != MPI_SUCCESS )
    throw std::runtime_error("MPI_Allreduce failed");
}

/**
 * @brief Function <code>AllReduce</code> copies the source/destination vector into a
 * temporary vector and then executed the MPI operation using the temporary as the source.
 *
 * @param mpi_comm	a <code>MPI_Comm</code> value of the MPI communicator.
 *
 * @param op		a <code>MPI_Op</code> value of the MPI operation.
 *
 * @param src_dest	a <code>std::vector<T></code> reference to be copied for the
 *			source and used as the destination.
 *
 * @param size		a <code>size_t</code> value of the length of the array pointed to
 *			by <b>src_dest</b>
 *
 */
template<class T>
inline void
AllReduce(MPI_Comm mpi_comm, MPI_Op op, std::vector<T> &dest)
{
  std::vector<T> source(dest);

  if (MPI_Allreduce(&source[0], &dest[0], (int) dest.size(), Datatype<T>::type(), op, mpi_comm) != MPI_SUCCESS )
    throw std::runtime_error("MPI_Allreduce failed");
}

/**
 * @brief Function <code>AllReduce</code> copies the source/destination vector into a
 * temporary vector and then executed the MPI operation using the temporary as the source.
 *
 * @param mpi_comm	a <code>MPI_Comm</code> value of the MPI communicator.
 *
 * @param op		a <code>MPI_Op</code> value of the MPI operation.
 *
 * @param src		a <code>std::vector<T></code> reference to the source data for the
 *			MPI op.
 *
 * @param dest		a <code>std::vector<T></code> reference to the destination data
 *			for the MPI op.
 *
 * @param size		a <code>size_t</code> value of the length of the array pointed to
 *			by <b>src_dest</b>
 *
 */
template<class T>
inline void
AllReduce(MPI_Comm mpi_comm, MPI_Op op, std::vector<T> &source, std::vector<T> &dest)
{
  if (source.size() != dest.size())
    throw std::runtime_error("sierra::MPI::AllReduce(MPI_Comm mpi_comm, MPI_Op op, std::vector<T> &source, std::vector<T> &dest) vector lengths not equal");

  if (MPI_Allreduce(&source[0], &dest[0], (int) dest.size(), Datatype<T>::type(), op, mpi_comm) != MPI_SUCCESS )
    throw std::runtime_error("MPI_Allreduce failed");
}


template<class T>
inline void
AllGather(MPI_Comm mpi_comm, std::vector<T> &source, std::vector<T> &dest)
{
  int nproc = 1;
  MPI_Comm_size(mpi_comm,&nproc);
  if (source.size()*nproc != dest.size())
    throw std::runtime_error("sierra::MPI::AllReduce(MPI_Comm mpi_comm, MPI_Op op, std::vector<T> &source, std::vector<T> &dest) vector lengths not equal");

  if (MPI_Allgather(&source[0], (int)source.size(), Datatype<T>::type(),
                    &dest[0],   (int)source.size(), Datatype<T>::type(),
                    mpi_comm) != MPI_SUCCESS ){
    throw std::runtime_error("MPI_Allreduce failed");
  }
}

/**
 * @brief Function <code>align_cast</code> returns a pointer that has been aligned to the
 * specified alignment or double if the alignment if greater than that of double.
 *
 * @param p		a <code>void</code> pointer of the address to align.
 *
 * @return		a <code>T</code> pointer which is an alignment of <b>p</b> to the
 *			lesser of the type specified or double.
 */
template <typename T>
T *align_cast(void *p)
{
  enum {alignment = (sizeof(T) > sizeof(double) ? sizeof(double) : sizeof(T))};
  enum {mask = alignment - 1};

  char * c = reinterpret_cast<char *>(p);
  size_t front_misalign = (c - (char *)0) & mask;
  if (front_misalign > 0) {
    size_t correction = alignment - front_misalign;
    T *q = reinterpret_cast<T *>((c - (char *)0) + correction);
    return q;
  }

  return reinterpret_cast<T *>(p);
}

/**
 * @brief Interface class <code>ReduceInterface</code> specifies the required virtual
 * functions for the aggregated type and operation operator.  The aggregated reduction
 * operator allows a single MPI operation to perform many operations on many types.
 *
 * This is accomplished by the creation of a user MPI op, which consists of a vector of
 * this interface.  This vector is used to first size the data, then to create source and
 * destination buffers, copy the source data into the source buffer, perform the operation
 * on each piece of data from the source to the destination buffer and finally copy the
 * destination buffer to destination data.
 *
 */
struct ReduceInterface {
  /**
   * Creates a new <b>ReduceInterface</b> instance.
   *
   */
  ReduceInterface()
  {}

  /**
   * Destroys a <b>ReduceInterface</b> instance.
   *
   */
  virtual ~ReduceInterface()
  {}

  /**
   * @brief Member function <b>size</b> returns the size in bytes needed to store the data for
   * the reduction operation of this interface.  The <b>inbuf</b> parameter should be
   * advanced to provide enough aligned space for the data of the reduction.
   *
   * @param inbuf		a <b>void</b> reference to a pointer to the be advanced by
   *				the size needed to store the data for the reduction.
   *
   */
  virtual void size(void *&inbuf) const = 0;

  /**
   * @brief Member function <b>copyin</b> copies the data from the reduction interface to
   * the <b>inbuf</b> pointer reference.  The pointer should be advanced to point just
   * beyond the space of the data.  Be sure to align the pointer before storing the first
   * element of data.
   *
   * @param inbuf		a <b>void</b> reference to a pointer to start placing data
   *				after alignment and should be advance to just beyond the
   *				last element of stored data.
   *
   */
  virtual void copyin(void *&inbuf) const = 0;

  /**
   * @brief Member function <b>copyin</b> copies the data from <b>outbuf</b> pointer
   * reference to the reduction interface data.  The pointer should be advanced to point
   * just beyond the space of the data.  Be sure to align the pointer before retrieving
   * the first element of data.
   *
   * @param outbuf		a <b>void</b> reference to a pointer to start retrieving
   *				data after alignment and should be advance to just beyond
   *				the last element of retrieved data.
   *
   */
  virtual void copyout(void *&outbuf) const = 0;

  /**
   * @brief Member function <b>op</b> executes the operation on the data at <b>inbuf</b>
   * pointer reference and <b>outbuf</b> pointer reference and placing the result in
   * <b>outbuf</b> pointer reference.  The pointers should be advanced to point just
   * beyond the space of the correspondingdata.  Be sure to align the pointer before
   * retrieving the first element of data.
   *
   * @param inbuf		a <b>void</b> reference to a pointer to start retrieving data
   *				after alignment and should be advance to just beyond the
   *				last element of stored data.
   *
   * @param outbuf		a <b>void</b> reference to a pointer to start retrieving
   *				and storing data after alignment and should be advance to
   *				just beyond the last element of retrieved data.
   *
   */
  virtual void op(void *&inbuf, void *&outbuf) const = 0;
};


/**
 * @brief Template class <code>Reduce</code> implements the ReduceInterface interface for
 * any operator and type.  The operator is a functor and the iterator specifies the access
 * to the data and its type.
 *
 * The operator <b>Op</b> is a functor that accepts two pointers, one to the destination
 * object of the iterator's value type and one to the source object of the iterator's
 * value type.  The <b>It</b> is an iterator for accessing the data.  Remember that a
 * pointer is an iterator, so using plain arrays, even of length one as in a pointer to a
 * scalar meets the criteria.
 *
 */
template<class Op, class LocalIt, class GlobalIt = LocalIt>
struct Reduce : public ReduceInterface
{
  typedef typename std::iterator_traits<LocalIt>::value_type value_type;
  typedef typename std::iterator_traits<LocalIt>::difference_type difference_type;

  Reduce(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end)
    : m_localBegin(local_begin),
      m_localEnd(local_end),
      m_globalBegin(global_begin),
      m_globalEnd(global_end),
      m_length(local_end - local_begin)
  {
    if (global_end - global_begin != m_length)
      throw std::runtime_error("sierra::MPI::Reduce::Reduce(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end) local and global lengths not equal");
  }

  virtual ~Reduce()
  {}

  virtual void size(void *&inbuf) const {
    value_type *t = align_cast<value_type>(inbuf);
    t += m_length;
    inbuf = t;
  }

  virtual void copyin(void *&inbuf) const {
    value_type *t = align_cast<value_type>(inbuf);
    for (LocalIt it = m_localBegin; it != m_localEnd; ++it)
      *t++ = (*it);
    inbuf = t;
  }

  virtual void copyout(void *&outbuf) const {
    value_type *t = align_cast<value_type>(outbuf);
    for (GlobalIt it = m_globalBegin; it != m_globalEnd; ++it)
      (*it) = *t++;
    outbuf = t;
  }

  virtual void op(void *&inbuf, void *&outbuf) const {
    value_type *tin = align_cast<value_type>(inbuf);
    value_type *tout = align_cast<value_type>(outbuf);

    for (size_t i = m_length; i; --i)
      Op(tout++, tin++);
    inbuf = tin;
    outbuf = tout;
  }

  LocalIt			m_localBegin;
  LocalIt			m_localEnd;
  GlobalIt			m_globalBegin;
  GlobalIt			m_globalEnd;
  difference_type		m_length;
};


/**
 * @brief Class <code>ReduceSet</code> ...
 *
 */
class ReduceSet
{
public:
  typedef std::vector<ReduceInterface *> ReduceVector;

  ReduceSet();

  virtual ~ReduceSet();

  void add(ReduceInterface *reduce_interface);

  size_t size() const;

  void copyin(void * const buffer_in) const;

  void copyout(void * const buffer_out) const;

  void op(void * const buffer_in, void * const buffer_out) const;

  static void void_op(void * inv, void * outv, int *n, MPI_Datatype *datatype);

private:
  ReduceVector		m_reduceVector;
};

/**
 * @brief Member function <code>AllReduce</code> ...
 *
 * @param comm		a <code>MPI_Comm</code> variable ...
 *
 * @param reduce_set	a <code>ReduceSet</code> variable ...
 *
 */
void AllReduce(MPI_Comm comm, const ReduceSet &reduce_set);

/**
 * @brief Class <code>Sum</code> ...
 *
 */
struct Sum
{
  template <typename T>
  inline Sum(T * dest, const T *source) {
    *dest += *source;
  }
};

/**
 * @brief Class <code>Prod</code> ...
 *
 */
struct Prod
{
  template <typename T>
  inline Prod(T * dest, const T *source) {
    *dest *= *source;
  }
};

/**
 * @brief Class <code>Min</code> ...
 *
 */
struct Min
{
  template <typename T>
  inline Min(T * dest, const T *source) {
    *dest = std::min(*dest, *source);
  }
};

/**
 * @brief Class <code>Max</code> ...
 *
 */
struct Max
{
  template <typename T>
  inline Max(T * dest, const T *source) {
    *dest = std::max(*dest, *source);
  }
};

/**
 * @brief Class <code>MinLoc</code> ...
 *
 */
struct MinLoc
{
  template <typename T>
  inline MinLoc(Loc<T> * dest, const Loc<T> *source) {
    if (source->m_value < dest->m_value) {
      dest->m_value = source->m_value;
      dest->m_loc = source->m_loc;
    }
    else if (source->m_value == dest->m_value)
      dest->m_loc = std::min(dest->m_loc, source->m_loc);
  }
};

/**
 * @brief Class <code>MaxLoc</code> ...
 *
 */
struct MaxLoc
{
  template <typename T>
  inline MaxLoc(Loc<T> * dest, const Loc<T> *source) {
    if (source->m_value > dest->m_value) {
      dest->m_value = source->m_value;
      dest->m_loc = source->m_loc;
    }
    else if (source->m_value == dest->m_value)
      dest->m_loc = std::min(dest->m_loc, source->m_loc);
  }
};


struct MaxTempLoc
{
  inline MaxTempLoc(TempLoc * dest, const TempLoc *source) {
    if (source->m_value > dest->m_value) {
      dest->m_value = source->m_value;
      dest->m_other = source->m_other;
      dest->m_loc = source->m_loc;
    }
    else if (source->m_value == dest->m_value) {
      if (dest->m_loc > source->m_loc) {
        dest->m_other = source->m_other;
        dest->m_loc = source->m_loc;
      }
    }
  }
};

struct MinTempLoc
{
  inline MinTempLoc(TempLoc * dest, const TempLoc *source) {
    if (source->m_value < dest->m_value) {
      dest->m_value = source->m_value;
      dest->m_other = source->m_other;
      dest->m_loc = source->m_loc;
    }
    else if (source->m_value == dest->m_value) {
      if (dest->m_loc > source->m_loc) {
        dest->m_other = source->m_other;
        dest->m_loc = source->m_loc;
      }
    }
  }
};

/**
 * @brief Member function <code>ReduceSum</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @param length	a <code>size_t</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Sum, T *> *ReduceSum(T *t, T *u, size_t length) {
  return new Reduce<Sum, T *>(t, t + length, u, u + length);
}

/**
 * @brief Member function <code>ReduceProd</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @param length	a <code>size_t</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Prod, T *> *ReduceProd(T *t, T *u, size_t length) {
  return new Reduce<Prod, T *>(t, t + length, u, u + length);
}

/**
 * @brief Member function <code>ReduceMax</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @param length	a <code>size_t</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Max, T *> *ReduceMax(T *t, T *u, size_t length) {
  return new Reduce<Max, T *>(t, t + length, u, u + length);
}

/**
 * @brief Member function <code>ReduceMin</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @param length	a <code>size_t</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Min, T *> *ReduceMin(T *t, T *u, size_t length) {
  return new Reduce<Min, T *>(t, t + length, u, u + length);
}


/**
 * @brief Member function <code>ReduceSum</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Sum, T *> *ReduceSum(T &t, T &u) {
  return new Reduce<Sum, T *>(&t, &t + 1, &u, &u + 1);
}

/**
 * @brief Member function <code>ReduceProd</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Prod, T *> *ReduceProd(T &t, T &u) {
  return new Reduce<Prod, T *>(&t, &t + 1, &u, &u + 1);
}

/**
 * @brief Member function <code>ReduceMax</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Max, T *> *ReduceMax(T &t, T &u) {
  return new Reduce<Max, T *>(&t, &t + 1, &u, &u + 1);
}

/**
 * @brief Member function <code>ReduceMin</code> ...
 *
 * @param t		a <code>T</code> variable ...
 *
 * @param u		a <code>T</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T>
Reduce<Min, T *> *ReduceMin(T &t, T &u) {
  return new Reduce<Min, T *>(&t, &t + 1, &u, &u + 1);
}


/**
 * @brief Member function <code>ReduceSum</code> ...
 *
 * @param local_begin	an <code>LocalIt</code> variable ...
 *
 * @param local_end	an <code>LocalIt</code> variable ...
 *
 * @param global_begin	an <code>GlobalIt</code> variable ...
 *
 * @param global_end	an <code>GlobalIt</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<class LocalIt, class GlobalIt>
Reduce<Sum, LocalIt, GlobalIt> *ReduceSum(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end) {
  return new Reduce<Sum, LocalIt, GlobalIt>(local_begin, local_end, global_begin, global_end);
}

/**
 * @brief Member function <code>ReduceProd</code> ...
 *
 * @param local_begin	an <code>LocalIt</code> variable ...
 *
 * @param local_end	an <code>LocalIt</code> variable ...
 *
 * @param global_begin	an <code>GlobalIt</code> variable ...
 *
 * @param global_end	an <code>GlobalIt</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<class LocalIt, class GlobalIt>
Reduce<Prod, LocalIt, GlobalIt> *ReduceProd(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end) {
  return new Reduce<Prod, LocalIt, GlobalIt>(local_begin, local_end, global_begin, global_end);
}

/**
 * @brief Member function <code>ReduceMin</code> ...
 *
 * @param local_begin	an <code>LocalIt</code> variable ...
 *
 * @param local_end	an <code>LocalIt</code> variable ...
 *
 * @param global_begin	an <code>GlobalIt</code> variable ...
 *
 * @param global_end	an <code>GlobalIt</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T, class LocalIt, class GlobalIt>
Reduce<Min, LocalIt, GlobalIt> *ReduceMin(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end) {
  return new Reduce<Min, LocalIt>(local_begin, local_end, global_begin, global_end);
}

/**
 * @brief Member function <code>ReduceMax</code> ...
 *
 * @param local_begin	an <code>LocalIt</code> variable ...
 *
 * @param local_end	an <code>LocalIt</code> variable ...
 *
 * @param global_begin	an <code>GlobalIt</code> variable ...
 *
 * @param global_end	an <code>GlobalIt</code> variable ...
 *
 * @return a <code>Reduce</code> ...
 */
template<typename T, class LocalIt, class GlobalIt>
Reduce<Max, LocalIt, GlobalIt> *ReduceMax(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end) {
  return new Reduce<Max, LocalIt>(local_begin, local_end, global_begin, global_end);
}

/**
 * @brief Member function <code>AllReduceCollected</code> ...
 *
 * @param mpi_comm	a <code>MPI_Comm</code> variable ...
 *
 * @param op		a <code>MPI_Op</code> variable ...
 *
 * @param collector	an <code>U</code> variable ...
 *
 */
template<class T, class U>
inline void
AllReduceCollected(MPI_Comm mpi_comm, MPI_Op op, U collector)
{
  std::vector<T> source;

  std::back_insert_iterator<std::vector<T> > source_inserter(source);
  collector.gather(op, source_inserter);

  int size = source.size();

#ifdef SIERRA_DEBUG
  //
  //  Check that the array lengths being reduces are all the same
  //
  int num_proc;
  int my_proc;
  MPI_Comm_size(mpi_comm, &num_proc);
  MPI_Comm_rank(mpi_comm, &my_proc);


  std::vector<int> local_array_len(num_proc, 0);
  local_array_len[my_proc] == size;
  std::vector<int> global_array_len(num_proc, 0);

  MPI_Allreduce(&local_array_len[0], &global_array_len[0], num_proc, MPI_INT, MPI_SUM, mpi_comm);

  for(unsigned i = 0; i < num_proc; ++i) {
    if(global_array_len[i] != size) {
      throw std::runtime_error("Slib_MPI.h::AllReduceCollected, not all processors have the same length array");
    }
  }
#endif

  if (source.empty()) return;
  std::vector<T> dest(size);

  if (MPI_Allreduce(&source[0], &dest[0], size, Datatype<T>::type(), op, mpi_comm) != MPI_SUCCESS )
    throw std::runtime_error("MPI_Allreduce failed");

  typename std::vector<T>::iterator dest_getter = dest.begin();
  collector.scatter(op, dest_getter);
}

///
/// @}
///

} // namespace MPI
} // namespace sierra

#endif // if defined( STK_HAS_MPI )
#endif // STK_UTIL_PARALLEL_MPI_hpp
