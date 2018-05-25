#include "mtk_kokkos.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Field.hpp"
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include "stk_mesh/base/FieldBLAS.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/getOption.h"

#include <sstream>

typedef stk::mesh::Field<double> ScalarField;

#ifdef KOKKOS_ENABLE_OPENMP
  typedef Kokkos::OpenMP   ExecSpace ;
#elif KOKKOS_ENABLE_CUDA
  typedef Kokkos::Cuda     ExecSpace ;
#else
  typedef Kokkos::Serial   ExecSpace ;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
   typedef Kokkos::OpenMP       MemSpace;
#elif KOKKOS_ENABLE_CUDA
   typedef Kokkos::CudaSpace    MemSpace;
#else
   typedef Kokkos::HostSpace    MemSpace;
#endif

  // typedef Kokkos::CudaUVMSpace MemSpace;

  typedef Kokkos::LayoutLeft   Layout ;
  // typedef Kokkos::LayoutRight  Layout ;

  typedef Kokkos::RangePolicy<ExecSpace> range_policy ;

  // Allocate y, x vectors and Matrix A:
  // Device
  typedef Kokkos::View<double*, Layout, MemSpace>   ViewVectorType;
  typedef ViewVectorType KokkosVector;

inline double elapsed_time(const struct timeval &begin, const struct timeval &end)
{
    return 1.0*(end.tv_sec-begin.tv_sec) + 1.0e-6*(end.tv_usec-begin.tv_usec);
}


template<class Scalar>
class NgpFieldVector
{
public:

    NgpFieldVector(stk::mesh::Field<Scalar>& field, const stk::mesh::Selector& selector) :
        scalarField(field),
        vector(field.name(), get_max_field_bucket_length(field, selector)),
        hostVector(Kokkos::create_mirror_view(vector))
    {
    }

    NgpFieldVector(stk::mesh::Field<Scalar>& field) :
        scalarField(field),
        vector(field.name(), get_max_field_bucket_length(field)),
        hostVector(Kokkos::create_mirror_view(vector))
    {
    }

    KokkosVector& device_get() {return vector;}
    KokkosVector::HostMirror& host_get() {return hostVector;}

    void copy_to_device(int kmax, Scalar x[])
    {
        // Initialize x vector on host
        for (int j = 0; j < kmax; ++j) {
          hostVector( j ) = x[j];
        }

        //Deep copy host view to device views
        Kokkos::deep_copy(vector, hostVector);
    }

    stk::mesh::Field<Scalar> &get_field() {return scalarField;}
    const stk::mesh::Field<Scalar> &get_field() const {return scalarField;}

    stk::mesh::BucketVector const& get_buckets()
    {
        return get_buckets(stk::mesh::selectField(scalarField));
    }
    stk::mesh::BucketVector const& get_buckets(const stk::mesh::Selector& selector)
    {
        return scalarField.get_mesh().get_buckets(scalarField.entity_rank(),
                                                  selector & scalarField.mesh_meta_data().locally_owned_part());
    }

private:

    unsigned get_max_field_bucket_length(const stk::mesh::Field<Scalar> & fieldBase)
    {
        return get_max_field_bucket_length(fieldBase, stk::mesh::selectField(fieldBase));
    }

    unsigned get_max_field_bucket_length(const stk::mesh::Field<Scalar> & fieldBase,
                                         const stk::mesh::Selector& selector)
    {
        stk::mesh::BucketVector const& buckets = fieldBase.get_mesh().get_buckets(fieldBase.entity_rank(),
                                                                                   selector & fieldBase.mesh_meta_data().locally_owned_part());

        unsigned max_field_bucket_length = 0;

        for(size_t i=0; i < buckets.size(); i++) {
            stk::mesh::Bucket & b = *buckets[i];
            const stk::mesh::Bucket::size_type length = b.size();
            const unsigned int fieldSize = field_scalars_per_entity(fieldBase, b);
            const unsigned kmax = length * fieldSize;
            max_field_bucket_length = std::max(max_field_bucket_length, kmax);
        }
        return max_field_bucket_length;
    }

private:
    stk::mesh::Field<Scalar> &scalarField;
    KokkosVector vector;
    KokkosVector::HostMirror hostVector;
};

template<class Scalar>
struct KokkosBLAS
{
    inline
    static void axpy( const int & kmax, const Scalar & alpha, const Scalar x[], Scalar y[])
    {
        for(int k = 0 ; k < kmax ; ++k) {
            y[k] = alpha * x[k] + y[k];
        }
    }

    inline
    static void copy( const int & kmax, const Scalar x[], Scalar y[])
    {
        for(int k = 0 ; k < kmax ; ++k) {
            y[k] = x[k];
        }
    }

    inline
    static void product( const int & kmax, const Scalar x[], const Scalar y[], Scalar z[])
    {
        for (int k=0;k<kmax;k++)
        {
            z[k]=x[k]*y[k];
        }
    }

    inline
    static Scalar dot( const int & kmax, const Scalar x[], const Scalar y[])
    {
        Scalar result = Scalar(0.0);
        for(int k = 0 ; k < kmax ; ++k) {
            result+=y[k] * x[k];
        }
        return result;
    }

    inline
    static Scalar nrm2( const int & kmax, const Scalar x[])
    {
        Scalar result = Scalar(0.0);
        for(int k = 0 ; k < kmax ; ++k) {
            result+=pow(std::abs(x[k]),2);
        }
        return Scalar(sqrt(result));
    }

    inline
    static void scal( const int & kmax, const Scalar alpha, Scalar x[])
    {
        for(int k = 0 ; k < kmax ; ++k) {
            x[k] = alpha * x[k];
        }
    }

    inline
    static void fill(const int & kmax, const Scalar alpha, Scalar x[],const int inc=1)
    {
        for(int k = 0 ; k < kmax*inc ; k+=inc) {
            x[k] = alpha;
        }
    }

    //    inline
    //    static void fill(const int & kmax, const Scalar alpha, Scalar x[])
    //    {
    //        std::fill(x,x+kmax,alpha);
    //    }

    inline
    static void swap(const int & kmax, Scalar x[], Scalar y[])
    {
        Scalar temp;
        for(int k = 0 ; k < kmax ; ++k) {
            temp = y[k];
            y[k] = x[k];
            x[k] = temp;
        }
    }

    inline
    static Scalar asum(const int & kmax, const Scalar x[])
    {
        Scalar result = Scalar(0.0);
        for(int k = 0 ; k < kmax ; ++k) {
            result+=std::abs(x[k]);
        }
        return Scalar(result);
    }

    inline
    static int iamax( const int & kmax, const Scalar x[])
    {
        double amax = 0.0;
        int result = 0;
        for(int k = 0 ; k < kmax ; ++k) {
            if (amax<std::abs(x[k])) {
                result = k;
                amax = std::abs(x[k]);
            }
        }
        return result;
    }

    inline
    static int iamin( const int & kmax, const Scalar x[])
    {
        int result = 0;
        double amin = std::abs(x[0]);
        for(int k = 0 ; k < kmax ; ++k) {
            if (std::abs(x[k])<amin) {
                result = k;
                amin = std::abs(x[k]);
            }
        }
        return result;
    }
};

namespace stk {
namespace kokkos {
    // **********************************************************************************************************************************************
    template<class Scalar>
    inline
    int field_length(NgpFieldVector<Scalar> & xField, const stk::mesh::Bucket & b)
    {
        const stk::mesh::Bucket::size_type length = b.size();
        const unsigned int fieldSize = field_scalars_per_entity(xField.get_field(), b);
        return length * fieldSize;
    }

    template<class Scalar>
    inline
    Scalar * field_ptr(NgpFieldVector<Scalar> & xField, const stk::mesh::Bucket & b)
    {
        return static_cast<Scalar*>(field_data(xField.get_field(), b));
    }
    // **********************************************************************************************************************************************
    template<class Scalar>
    inline
    void field_fill(
            const Scalar data,
            NgpFieldVector<Scalar> & xField,
            const stk::mesh::Selector selector,
            const MPI_Comm comm)
    {
        ThrowAssert(typeid(Scalar) == xField.get_field().data_traits().type_info);

        stk::mesh::BucketVector const& buckets = xField.get_buckets(selector);

        for(size_t i=0; i < buckets.size(); i++) {
            stk::mesh::Bucket & b = *buckets[i];
            KokkosBLAS<Scalar>::fill(field_length(xField, b), data, field_ptr(xField, b));
        }
    }

    template<class Scalar>
    inline
    void field_fill(
            const Scalar data,
            NgpFieldVector<Scalar> & xField,
            const stk::mesh::Selector selector)
    {
        const MPI_Comm comm = xField.get_field().get_mesh().parallel();
        stk::kokkos::field_fill(data,xField,selector,comm);
    }

    template<class Scalar>
    inline
    void field_fill(
            const Scalar data,
            NgpFieldVector<Scalar> & xField)
    {
        const stk::mesh::Selector selector = stk::mesh::selectField(xField.get_field());
        stk::kokkos::field_fill(data,xField,selector);
    }
    // **********************************************************************************************************************************************
    template<class Scalar>
    inline
    void field_asum(
            Scalar & glob_result,
            NgpFieldVector<Scalar>& xField,
            const stk::mesh::Selector& selector,
            const MPI_Comm comm)
    {
        ThrowAssert(typeid(Scalar) == xField.get_field().data_traits().type_info);

        stk::mesh::BucketVector const& buckets = xField.get_buckets(selector);
        KokkosVector& dev_x = xField.device_get();

        Scalar local_result = Scalar(0.0);

        for(size_t i=0; i < buckets.size(); i++) {
            stk::mesh::Bucket & b = *buckets[i];
            int kmax = field_length(xField, b);
            xField.copy_to_device(kmax, field_ptr(xField, b));
            Scalar gpu_result = Scalar(0.0);

            Kokkos::parallel_reduce( range_policy( 0, kmax), KOKKOS_LAMBDA ( int k, double& update ) {
              update += dev_x( k );
            }, gpu_result );

            local_result += gpu_result;
        }

        glob_result=local_result;
        stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    }

    template<class Scalar>
    inline
    void field_asum(
            Scalar & result,
            NgpFieldVector<Scalar>& xFieldBase,
            const stk::mesh::Selector& selector)
    {
        const MPI_Comm comm = xFieldBase.get_field().get_mesh().parallel();
        stk::kokkos::field_asum(result,xFieldBase,selector,comm);
    }

    template<class Scalar>
    inline
    void field_asum(
            Scalar & result,
            NgpFieldVector<Scalar>& xFieldBase)
    {
        const stk::mesh::Selector selector = stk::mesh::selectField(xFieldBase.get_field());
        stk::kokkos::field_asum(result,xFieldBase,selector);
    }
    // **********************************************************************************************************************************************
    template<class Scalar>
    inline
    void field_nrm2(
            Scalar & glob_result,
            NgpFieldVector<Scalar>& xField,
            const stk::mesh::Selector& selector,
            const MPI_Comm comm)
    {
        ThrowAssert(typeid(Scalar) == xField.get_field().data_traits().type_info);

        stk::mesh::BucketVector const& buckets = xField.get_buckets(selector);
        KokkosVector& dev_x = xField.device_get();

        Scalar local_result = Scalar(0.0);

        for(size_t i=0; i < buckets.size(); i++) {
            stk::mesh::Bucket & b = *buckets[i];
            int kmax = field_length(xField, b);
            xField.copy_to_device(kmax, field_ptr(xField, b));
            Scalar gpu_result = Scalar(0.0);

            Kokkos::parallel_reduce( range_policy( 0, kmax), KOKKOS_LAMBDA ( int k, double& update ) {
              update += dev_x( k ) * dev_x( k );
            }, gpu_result );

            local_result += gpu_result;
        }

        glob_result=local_result;
        stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
        glob_result=sqrt(glob_result);
    }
    template<class Scalar>
    inline
    void field_nrm2(
            Scalar & result,
            NgpFieldVector<Scalar>& xField,
            const stk::mesh::Selector& selector)
    {
        const MPI_Comm comm = xField.get_field().get_mesh().parallel();
        stk::kokkos::field_nrm2(result,xField,selector,comm);
    }

    template<class Scalar>
    inline
    void field_nrm2(
            Scalar & result,
            NgpFieldVector<Scalar>& xField)
    {
        const stk::mesh::Selector selector = stk::mesh::selectField(xField.get_field());
        stk::kokkos::field_nrm2(result,xField,selector);
    }
    // **********************************************************************************************************************************************

    template<class Scalar>
    inline
    void field_dot(
            Scalar & glob_result,
            NgpFieldVector<Scalar>& xField,
            NgpFieldVector<Scalar>& yField,
            const stk::mesh::Selector& selector,
            const MPI_Comm comm)
    {
        ThrowAssert(xField.get_field().entity_rank() == yField.get_field().entity_rank());
        ThrowAssert(&xField.get_field().get_mesh() == &yField.get_field().get_mesh());
        ThrowAssert(xField.get_field().data_traits().type_info == yField.get_field().data_traits().type_info);
        ThrowAssert(typeid(Scalar) == xField.get_field().data_traits().type_info);

        stk::mesh::BucketVector const& buckets = xField.get_buckets(selector);
        KokkosVector& dev_x = xField.device_get();
        KokkosVector& dev_y = yField.device_get();

        Scalar local_result = Scalar(0.0);

        for(size_t i=0; i < buckets.size(); i++) {
            stk::mesh::Bucket & b = *buckets[i];
            int kmax = field_length(xField, b);
            xField.copy_to_device(kmax, field_ptr(xField, b));
            yField.copy_to_device(kmax, field_ptr(yField, b));
            Scalar gpu_result = Scalar(0.0);

            Kokkos::parallel_reduce( range_policy( 0, kmax), KOKKOS_LAMBDA ( int k, double& update ) {
              update += dev_x( k ) * dev_y( k );
            }, gpu_result );

            local_result += gpu_result;
        }

        glob_result = local_result;
        stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    }

    template<class Scalar>
    inline
    void field_dot(
            Scalar & result,
            NgpFieldVector<Scalar>& xField,
            NgpFieldVector<Scalar>& yField,
            const stk::mesh::Selector& selector)
    {
        const MPI_Comm comm = xField.get_field().get_mesh().parallel();
        stk::kokkos::field_dot(result,xField,yField,selector,comm);
    }

    template<class Scalar>
    inline
    void field_dot(
            Scalar & result,
            NgpFieldVector<Scalar>& xField,
            NgpFieldVector<Scalar>& yField)
    {
        const stk::mesh::Selector selector = stk::mesh::selectField(xField.get_field()) & stk::mesh::selectField(yField.get_field());
        stk::kokkos::field_dot(result,xField,yField,selector);
    }
    // **********************************************************************************************************************************************
}
}

void get_mesh_dimensions(int &nx, int &ny, int &nz)
{
    nx = stk::unit_test_util::get_command_line_option<int>("-nx", 10);
    ny = stk::unit_test_util::get_command_line_option<int>("-ny", 10);
    nz = stk::unit_test_util::get_command_line_option<int>("-nz", 10);
}

void populate_mesh(stk::mesh::BulkData& bulk, int nx, int ny, int nz)
{
    stk::mesh::MetaData &meta = bulk.mesh_meta_data();

    std::ostringstream os;
    os << "generated:" << nx << "x" << ny << "x" << nz;

    double val = 0.0;
    ScalarField *x = &meta.declare_field<ScalarField>(stk::topology::ELEMENT_RANK, "x");
    ScalarField *y = &meta.declare_field<ScalarField>(stk::topology::ELEMENT_RANK, "y");

    stk::mesh::put_field(*x, meta.locally_owned_part(), &val);
    stk::mesh::put_field(*y, meta.locally_owned_part(), &val);

    stk::io::fill_mesh(os.str(), bulk);
}

TEST_F(MTK_Kokkos, stkFieldBLAS) {

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA);

    int nx, ny, nz;
    get_mesh_dimensions(nx, ny, nz);
    int dim = nx*ny*nz;

    populate_mesh(bulk, nx, ny, nz);

    struct timeval begin,end;
    double time, result, xResult, yResult;
    double Gbytes = double(sizeof(double) * ( 2* dim ))/(1024.0 * 1024.0 * 1024.0);

    ScalarField *x = static_cast<ScalarField *>(meta.get_field(stk::topology::ELEMENT_RANK, "x"));
    ScalarField *y = static_cast<ScalarField *>(meta.get_field(stk::topology::ELEMENT_RANK, "y"));

    {
        gettimeofday(&begin,NULL);
        double alpha = 1.0;
        stk::mesh::field_fill(alpha, *x);
        stk::mesh::field_fill(alpha, *y);
        gettimeofday(&end,NULL);

        std::cerr << "fill Time: " << elapsed_time(begin, end) << std::endl;
    }

    {
        gettimeofday(&begin,NULL);
        stk::mesh::field_asum(xResult, *x);
        stk::mesh::field_asum(yResult, *y);
        gettimeofday(&end,NULL);
        time = elapsed_time(begin, end);
        std::cerr << "asum Time: " << time << " bandwidth( " << Gbytes/time << " GB/s )" << std::endl;
        EXPECT_EQ(xResult, dim);
        EXPECT_EQ(yResult, dim);
    }

    {
        gettimeofday(&begin,NULL);
        stk::mesh::field_nrm2(xResult, *x);
        stk::mesh::field_nrm2(yResult, *y);
        gettimeofday(&end,NULL);
        time = elapsed_time(begin, end);
        std::cerr << "nrm2 Time: " << time << " bandwidth( " << Gbytes/time << " GB/s )" << std::endl;
        EXPECT_EQ(xResult, std::sqrt(1.0*dim));
        EXPECT_EQ(yResult, std::sqrt(1.0*dim));
    }

    {
        gettimeofday(&begin,NULL);
        stk::mesh::field_dot(result, *x, *y);
        gettimeofday(&end,NULL);
        time = elapsed_time(begin, end);
        std::cerr << "dot Time: " << time << " bandwidth( " << Gbytes/time << " GB/s )" << std::endl;
        EXPECT_EQ(result, dim);
    }
}



TEST_F(MTK_Kokkos, kokkosFieldBLAS) {

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA);

    int nx, ny, nz;
    get_mesh_dimensions(nx, ny, nz);
    int dim = nx*ny*nz;

    populate_mesh(bulk, nx, ny, nz);

    ScalarField *xField = static_cast<ScalarField *>(meta.get_field(stk::topology::ELEMENT_RANK, "x"));
    ScalarField *yField = static_cast<ScalarField *>(meta.get_field(stk::topology::ELEMENT_RANK, "y"));

    // This must be initialized after field is populated
    NgpFieldVector<double> x(*xField);
    NgpFieldVector<double> y(*yField);

    struct timeval begin,end;
    double time, result, xResult, yResult;
    double Gbytes = double(sizeof(double) * ( 2* dim ))/(1024.0 * 1024.0 * 1024.0);

    {
        gettimeofday(&begin,NULL);
        double alpha = 1.0;
        stk::kokkos::field_fill(alpha, x);
        stk::kokkos::field_fill(alpha, y);
        gettimeofday(&end,NULL);

        std::cerr << "fill Time: " << elapsed_time(begin, end) << std::endl;
    }

    {
        gettimeofday(&begin,NULL);
        stk::kokkos::field_asum(xResult, x);
        stk::kokkos::field_asum(yResult, y);
        gettimeofday(&end,NULL);
        time = elapsed_time(begin, end);
        std::cerr << "asum Time: " << time << " bandwidth( " << Gbytes/time << " GB/s )" << std::endl;
        EXPECT_EQ(xResult, dim);
        EXPECT_EQ(yResult, dim);
    }

    {
        gettimeofday(&begin,NULL);
        stk::kokkos::field_nrm2(xResult, x);
        stk::kokkos::field_nrm2(yResult, y);
        gettimeofday(&end,NULL);
        time = elapsed_time(begin, end);
        std::cerr << "nrm2 Time: " << time << " bandwidth( " << Gbytes/time << " GB/s )" << std::endl;
        EXPECT_EQ(xResult, std::sqrt(1.0*dim));
        EXPECT_EQ(yResult, std::sqrt(1.0*dim));
    }

    {
        gettimeofday(&begin,NULL);
        stk::kokkos::field_dot(result, x, y);
        gettimeofday(&end,NULL);
        time = elapsed_time(begin, end);
        std::cerr << "dot Time: " << time << " bandwidth( " << Gbytes/time << " GB/s )" << std::endl;
        EXPECT_EQ(result, dim);
    }
}
