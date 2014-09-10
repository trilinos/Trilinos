/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <stk_mesh/base/FieldBLAS.hpp>

#include <complex>
#include <string>

#if defined(_OPENMP) && !defined(__INTEL_COMPILER)
#define OPEN_MP_ACTIVE_PERFOPENMP_CPP
// there seems to be an issue with OpenMP combined with GoogleTest macros with Intel compilers
// example of error:
//    openMP.C(206): internal error: assertion failed at: "shared/cfe/edgcpfe/checkdir.c", line 5531
#include <omp.h>
#endif

//#ifndef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//TEST(perf_testing_openmp,failing) {
//    std::cout << "!defined(_OPENMP) || !defined(__INTEL_COMPILER)" <<std::endl;
//    EXPECT_TRUE(false);
//}
//#endif
//
//template<class Scalar_in,class Scalar_out>
//class TestFixture {
//private:
//    double test_time_0;
//    double test_time_1;
//
//    Scalar_out (*func_0_0)(const Scalar_in*,const int);
//    Scalar_out (*func_0_1)(const Scalar_in*,const unsigned int,const unsigned int);
//
//    Scalar_out (*func_1_0)(const Scalar_in*,const int);
//    Scalar_out (*func_1_1)(const Scalar_in*,const unsigned int,const unsigned int);
//public:
//    double get_time(unsigned int);
//    void return_to_the_beginning_of_time();
//
//    void set_func_0_0(Scalar_out (*func_0_tmp)(const Scalar_in*,const int));
//    void set_func_0_1(Scalar_out (*func_1_tmp)(const Scalar_in*,const unsigned int,const unsigned int));
//
//    void set_func_1_0(Scalar_out (*func_0_tmp)(const Scalar_in*,const Scalar_in*,const int));
//    void set_func_1_1(Scalar_out (*func_1_tmp)(const Scalar_in*,const Scalar_in*,const unsigned int,const unsigned int));
//
//    Scalar_out call_func_0(const Scalar_in*,const int);
//    Scalar_out call_func_1(const Scalar_in*,const unsigned int,const unsigned int);
//
//    std::string func_name;
//
//    TestFixture()
//    :test_time_0(0.0),
//     test_time_1(0.0) {}
//};
//
//template<class Scalar_in,class Scalar_out>
//void TestFixture<Scalar_in,Scalar_out>::return_to_the_beginning_of_time()
//{
//    test_time_0=0.0;
//    test_time_1=0.0;
//}
//
//template<class Scalar_in,class Scalar_out>
//double TestFixture<Scalar_in,Scalar_out>::get_time(unsigned int i)
//{
//    if (i==0) { return test_time_0; }
//    else if (i==1u) { return test_time_1; }
//    return -1.0;
//}
//
//template<class Scalar_in,class Scalar_out>
//void TestFixture<Scalar_in,Scalar_out>::set_func_1_0(Scalar_out (*func_0_tmp)(const Scalar_in*,const Scalar_in*,const int))
//{
//    func_1_0=func_0_tmp;
//}
//
//template<class Scalar_in,class Scalar_out>
//void TestFixture<Scalar_in,Scalar_out>::set_func_1_1(Scalar_out (*func_1_tmp)(const Scalar_in*,const Scalar_in*,const unsigned int,const unsigned int))
//{
//    func_1_1=func_1_tmp;
//}
//
//template<class Scalar_in,class Scalar_out>
//void TestFixture<Scalar_in,Scalar_out>::set_func_0_0(Scalar_out (*func_0_tmp)(const Scalar_in*,const int))
//{
//    func_0_0=func_0_tmp;
//}
//
//template<class Scalar_in,class Scalar_out>
//void TestFixture<Scalar_in,Scalar_out>::set_func_0_1(Scalar_out (*func_1_tmp)(const Scalar_in*,const unsigned int,const unsigned int))
//{
//    func_0_1=func_1_tmp;
//}
//
//template<class Scalar_in,class Scalar_out>
//Scalar_out TestFixture<Scalar_in,Scalar_out>::call_func_0(const Scalar_in* data,const int data_size)
//{
//    test_time_0-=stk::wall_time();
//    double result_tmp = func_0_0(data,data_size);
//    test_time_0+=stk::wall_time();
//    return result_tmp;
//}
//
//template<class Scalar_in,class Scalar_out>
//Scalar_out TestFixture<Scalar_in,Scalar_out>::call_func_1(const Scalar_in* data,const unsigned int num_buckets,const unsigned int bucket_size)
//{
//    test_time_1-=stk::wall_time();
//    double result_tmp = func_0_1(data,num_buckets,bucket_size);
//    test_time_1+=stk::wall_time();
//    return result_tmp;
//}
//
//double my_sum_0(const double* data,const int size)
//{
//    double result=0.0;
//    for (int i=0;i<size;i++)
//    {
//        result+=data[i];
//    }
//    return result;
//}
//
//double my_sum_1(const double* data,const unsigned int num_buckets,const unsigned int bucket_size)
//{
//    double result=0.0;
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        result+=my_sum_0(&data[i*bucket_size],bucket_size);
//    }
//    return result;
//}
//
//TEST(TestFixture,base)
//{
//    TestFixture<double,double> my_test_fixture;
//    std::string my_name=std::string("my func name");
//    my_test_fixture.func_name=my_name;
//
//    const double internal_value=5.21;
//
//    const unsigned int num_buckets=50;
//    const unsigned int bucket_size=512;
//    const int data_size=int(num_buckets*bucket_size);
//
//    double data [data_size];
//    for (unsigned int i=0;i<(unsigned int)data_size;i++)
//    {
//        data[i]=internal_value;
//    }
//
//    my_test_fixture.set_func_0_0(&my_sum_0);
//    my_test_fixture.set_func_0_1(&my_sum_1);
//
//    EXPECT_EQ(my_test_fixture.get_time(0),0.0);
//    EXPECT_EQ(my_test_fixture.get_time(1),0.0);
//
//    EXPECT_NEAR(data_size*internal_value,my_test_fixture.call_func_0(data,data_size),1.0e-5);
//
//    EXPECT_GT(my_test_fixture.get_time(0),0.0);
//    EXPECT_EQ(my_test_fixture.get_time(1),0.0);
//
//    EXPECT_NEAR(data_size*internal_value,my_test_fixture.call_func_1(data,num_buckets,bucket_size),1.0e-5);
//
//    EXPECT_GT(my_test_fixture.get_time(0),0.0);
//    EXPECT_GT(my_test_fixture.get_time(1),0.0);
//
//    my_test_fixture.return_to_the_beginning_of_time();
//
//    EXPECT_EQ(my_test_fixture.get_time(0),0.0);
//    EXPECT_EQ(my_test_fixture.get_time(1),0.0);
//
//    EXPECT_TRUE(my_name==my_test_fixture.func_name);
//}
//
//void perfOpenMP_general_test(void (*test_func)(unsigned int))
//{
//    const unsigned int len_num_buckets_list = 4u;
//    unsigned int num_buckets_list [len_num_buckets_list] = {4000u,2500u,1000u,500u};//{500u,1000u,2500u,4000u};
//    for (unsigned int k=0;k<len_num_buckets_list;k++)
//    {
//        test_func(num_buckets_list[k]);
//    }
//}
//
//double min_testing_serial(const double data[], const int data_size)
//{
//    double min = std::abs(data[0]);
//    double tmp;
//    for (int i=1;i<data_size;i++)
//    {
//        tmp=std::abs(data[i]);
//        if (tmp<min) min=tmp;
//    }
//    return min;
//}
//
//double min_testing_hiddenCritical(const double data[], const int data_size)
//{
//    double shar_min=std::abs(data[0]);
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) private(priv_tmp)
//#endif
//    for (int i=1;i<data_size;i++)
//    {
//        priv_tmp=std::abs(data[i]);
//        if (priv_tmp<shar_min)
//        {
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//            if (priv_tmp<shar_min) shar_min=priv_tmp;
//        }
//    }
//    return shar_min;
//}
//
//double min_testing_comparePrivates(const double data[], const int data_size)
//{
//    double shar_min=std::abs(data[0]);
//    double priv_min;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_min,priv_tmp)
//    {
//#endif
//        priv_min=std::abs(data[0]);
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (int i=1;i<data_size;i++)
//        {
//            priv_tmp=std::abs(data[i]);
//            if (priv_tmp<priv_min) priv_min=priv_tmp;
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (priv_min<shar_min)
//        {
//#pragma omp critical
//            if (priv_min<shar_min)
//#endif
//                shar_min=priv_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        }
//    }
//#endif
//    return shar_min;
//}
//
//double min_testing_reduction(const double data[], const int data_size)
//{
//    double shar_min=std::abs(data[0]);
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) reduction(min:shar_min) private(priv_tmp)
//#endif
//    for (int i=0;i<data_size;i++)
//    {
//        priv_tmp=std::abs(data[i]);
//        if (priv_tmp<shar_min) shar_min=priv_tmp;
//    }
//    return shar_min;
//}
//
//double min_testing_serial_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double min=std::abs(data[0]);
//    double tmp_min;
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        tmp_min=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size])]);
//        if (tmp_min<min) min=tmp_min;
//    }
//    return min;
//}
//
//double min_testing_hiddenCritical_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_min=std::abs(data[0]);
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) private(priv_tmp)
//#endif
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        priv_tmp=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size])]);
//        if (priv_tmp<shar_min)
//        {
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//            if (priv_tmp<shar_min) shar_min=priv_tmp;
//        }
//    }
//    return shar_min;
//}
//
//double min_testing_comparePrivates_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_min=std::abs(data[0]);
//    double priv_tmp;
//    double priv_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_tmp,priv_min)
//    {
//#endif
//        priv_min=std::abs(data[0]);
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (unsigned int i=0;i<num_buckets;i++)
//        {
//            priv_tmp=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size])]);
//            if (priv_tmp<priv_min)
//            {
//                priv_min=priv_tmp;
//            }
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (priv_min<shar_min)
//        {
//#pragma omp critical
//            if (priv_min<shar_min)
//#endif
//                shar_min=priv_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        }
//    }
//#endif
//    return shar_min;
//}
//
//double min_testing_reduction_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_min=std::abs(data[0]);
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) reduction(min:shar_min) private(priv_tmp)
//#endif
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        priv_tmp=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size])]);
//        if (priv_tmp<shar_min) shar_min=priv_tmp;
//    }
//    return shar_min;
//}
//
//void perf_testing_openmp_min(const unsigned int num_buckets)
//{
//    std::cout<<std::endl<<"num_buckets,"<<num_buckets<<std::endl;
//    const double high_value = -4.27;
//    const double low_value = -3.14;
//    const unsigned int bucket_size = 512;
//    const int data_size = bucket_size*num_buckets;
//
//    const int num_numThreads = 4;
//    const int numThread_list [num_numThreads] = {1,10,20,30};
//
//    const unsigned int test_repeat = 750;
//
//    double * data = new double [data_size];
//    for (unsigned int i=0;i<(unsigned int)data_size;i++)
//    {
//        data[i]=high_value;
//    }
//
//    const unsigned int seperation_increment = 100;
//
//    for (unsigned int i=seperation_increment;i<(unsigned int)(data_size-data_size%(int)seperation_increment);i+=seperation_increment)
//    {
//        data[i]=high_value+(low_value-high_value)*double(int(i)-(int)seperation_increment)/
//                double(-2*(int)seperation_increment+data_size-data_size%seperation_increment);
//    }
//
//    unsigned int num_test_types = 4u;
//    TestFixture<double,double> test_fixture[num_test_types];
//    test_fixture[0].set_func_0_0(&min_testing_serial);
//    test_fixture[0].set_func_0_1(&min_testing_serial_buckets);
//    test_fixture[0].func_name=std::string ("serial");
//
//    test_fixture[1].set_func_0_0(&min_testing_hiddenCritical);
//    test_fixture[1].set_func_0_1(&min_testing_hiddenCritical_buckets);
//    test_fixture[1].func_name=std::string ("hiddenCritical");
//
//    test_fixture[2].set_func_0_0(&min_testing_comparePrivates);
//    test_fixture[2].set_func_0_1(&min_testing_comparePrivates_buckets);
//    test_fixture[2].func_name=std::string ("comparePrivates");
//
//    test_fixture[3].set_func_0_0(&min_testing_reduction);
//    test_fixture[3].set_func_0_1(&min_testing_reduction_buckets);
//    test_fixture[3].func_name=std::string ("reduction");
//
//    char coutBuffer [100];
//    sprintf(coutBuffer,"%15s , %10s , %10s , %10s , %10s , %10s\n","algorithm name","numThreads","NonBucket","NonBucket\%","Buckets","Buckets\%");
//    std::cout<<coutBuffer;
//
//    for (unsigned int j=0;j<num_test_types;j++)
//    {
//        double single_threaded_time_0 = 50.0;
//        double single_threaded_time_1 = 50.0;
//
//        for (int t=0;t<num_numThreads;t++)
//        {
//            if (j==0 && t>0) continue;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            omp_set_num_threads(numThread_list[t]);
//#endif
//            for (unsigned int tr=0;tr<test_repeat;tr++)
//            {
//                EXPECT_EQ(test_fixture[j].call_func_0(data,data_size),std::abs(low_value));
//                EXPECT_EQ(test_fixture[j].call_func_1(data,num_buckets,bucket_size),std::abs(low_value));
//            }
//
//            if (t==0)
//            {
//                single_threaded_time_0=test_fixture[j].get_time(0);
//                single_threaded_time_1=test_fixture[j].get_time(1u);
//            }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            sprintf(coutBuffer,"%15s , %10d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#else
//            sprintf(coutBuffer,"%15s , %3d of %3d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    0,numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#endif
//            std::cout<<coutBuffer;
//
//            test_fixture[j].return_to_the_beginning_of_time();
//        }
//    }
//
//    delete data;
//}
//
//TEST(perf_testing_openmp,min)
//{
//    perfOpenMP_general_test(&perf_testing_openmp_min);
//}
//
//double max_testing_serial(const double data[], const int data_size)
//{
//    double max = 0.0;
//    double tmp;
//    for (int i=1;i<data_size;i++)
//    {
//        tmp=std::abs(data[i]);
//        if (max<tmp) max=tmp;
//    }
//    return max;
//}
//
//double max_testing_hiddenCritical(const double data[], const int data_size)
//{
//    double shar_max=0.0;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) private(priv_tmp)
//#endif
//    for (int i=1;i<data_size;i++)
//    {
//        priv_tmp=std::abs(data[i]);
//        if (shar_max<priv_tmp)
//        {
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//            if (shar_max<priv_tmp) shar_max=priv_tmp;
//        }
//    }
//    return shar_max;
//}
//
//double max_testing_comparePrivates(const double data[], const int data_size)
//{
//    double shar_max=0.0;
//    double priv_max;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_max,priv_tmp)
//    {
//#endif
//        priv_max=0.0;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (int i=1;i<data_size;i++)
//        {
//            priv_tmp=std::abs(data[i]);
//            if (priv_max<priv_tmp) priv_max=priv_tmp;
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (shar_max<priv_max)
//        {
//#pragma omp critical
//            if (shar_max<priv_max)
//#endif
//                shar_max=priv_max;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        }
//    }
//#endif
//    return shar_max;
//}
//
//double max_testing_reduction(const double data[], const int data_size)
//{
//    double shar_max=0.0;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) reduction(max:shar_max) private(priv_tmp)
//#endif
//    for (int i=0;i<data_size;i++)
//    {
//        priv_tmp=std::abs(data[i]);
//        if (shar_max<priv_tmp) shar_max=priv_tmp;
//    }
//    return shar_max;
//}
//
//double max_testing_serial_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double max=0.0;
//    double tmp_max;
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        tmp_max=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamax((int)bucket_size,&data[i*bucket_size])]);
//        if (max<tmp_max) max=tmp_max;
//    }
//    return max;
//}
//
//double max_testing_hiddenCritical_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_max=0.0;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) private(priv_tmp)
//#endif
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        priv_tmp=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamax((int)bucket_size,&data[i*bucket_size])]);
//        if (shar_max<priv_tmp)
//        {
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//            if (shar_max<priv_tmp) shar_max=priv_tmp;
//        }
//    }
//    return shar_max;
//}
//
//double max_testing_comparePrivates_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_max=0.0;
//    double priv_tmp;
//    double priv_max;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_tmp,priv_max)
//    {
//#endif
//        priv_max=0.0;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (unsigned int i=0;i<num_buckets;i++)
//        {
//            priv_tmp=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamax((int)bucket_size,&data[i*bucket_size])]);
//            if (priv_max<priv_tmp)
//            {
//                priv_max=priv_tmp;
//            }
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (shar_max<priv_max)
//        {
//#pragma omp critical
//            if (shar_max<priv_max)
//#endif
//                shar_max=priv_max;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        }
//    }
//#endif
//    return shar_max;
//}
//
//double max_testing_reduction_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_max=0.0;
//    double tmp_priv;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) reduction(max:shar_max) private(tmp_priv)
//#endif
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        tmp_priv=std::abs(data[int(i*bucket_size)+stk::mesh::FortranBLAS<double>::iamax((int)bucket_size,&data[i*bucket_size])]);
//        if (shar_max<tmp_priv) shar_max=tmp_priv;
//    }
//    return shar_max;
//}
//
//void perf_testing_openmp_max(const unsigned int num_buckets)
//{
//    std::cout<<std::endl<<"num_buckets,"<<num_buckets<<std::endl;
//    const double low_value = -2.07;
//    const double high_value = -4.27;
//    const unsigned int bucket_size = 512;
//    const int data_size = bucket_size*num_buckets;
//
//    const int num_numThreads = 4;
//    const int numThread_list [num_numThreads] = {1,10,20,30};
//
//    const unsigned int test_repeat = 750;
//
//    double * data = new double [data_size];
//    for (unsigned int i=0;i<(unsigned int)data_size;i++)
//    {
//        data[i]=low_value;
//    }
//
//    const unsigned int seperation_increment = 100;
//
//    for (unsigned int i=seperation_increment;i<(unsigned int)(data_size-data_size%(int)seperation_increment);i+=seperation_increment)
//    {
//        data[i]=low_value+(high_value-low_value)*double(int(i)-(int)seperation_increment)/
//                double(-2*(int)seperation_increment+data_size-data_size%seperation_increment);
//    }
//
//    unsigned int num_test_types = 4u;
//    TestFixture<double,double> test_fixture[num_test_types];
//    test_fixture[0].set_func_0_0(&max_testing_serial);
//    test_fixture[0].set_func_0_1(&max_testing_serial_buckets);
//    test_fixture[0].func_name=std::string ("serial");
//
//    test_fixture[1].set_func_0_0(&max_testing_hiddenCritical);
//    test_fixture[1].set_func_0_1(&max_testing_hiddenCritical_buckets);
//    test_fixture[1].func_name=std::string ("hiddenCritical");
//
//    test_fixture[2].set_func_0_0(&max_testing_comparePrivates);
//    test_fixture[2].set_func_0_1(&max_testing_comparePrivates_buckets);
//    test_fixture[2].func_name=std::string ("comparePrivates");
//
//    test_fixture[3].set_func_0_0(&max_testing_reduction);
//    test_fixture[3].set_func_0_1(&max_testing_reduction_buckets);
//    test_fixture[3].func_name=std::string ("reduction");
//
//    char coutBuffer [100];
//    sprintf(coutBuffer,"%15s , %10s , %10s , %10s , %10s , %10s\n","algorithm name","numThreads","NonBucket","NonBucket\%","Buckets","Buckets\%");
//    std::cout<<coutBuffer;
//
//    for (unsigned int j=0;j<num_test_types;j++)
//    {
//        double single_threaded_time_0 = 50.0;
//        double single_threaded_time_1 = 50.0;
//
//        for (int t=0;t<num_numThreads;t++)
//        {
//            if (j==0 && t>0) continue;
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            omp_set_num_threads(numThread_list[t]);
//#endif
//            for (unsigned int tr=0;tr<test_repeat;tr++)
//            {
//                EXPECT_EQ(test_fixture[j].call_func_0(data,data_size),std::abs(high_value));
//                EXPECT_EQ(test_fixture[j].call_func_1(data,num_buckets,bucket_size),std::abs(high_value));
//            }
//
//            if (t==0)
//            {
//                single_threaded_time_0=test_fixture[j].get_time(0);
//                single_threaded_time_1=test_fixture[j].get_time(1u);
//            }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            sprintf(coutBuffer,"%15s , %10d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#else
//            sprintf(coutBuffer,"%15s , %3d of %3d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    0,numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#endif
//
//            std::cout<<coutBuffer;
//
//            test_fixture[j].return_to_the_beginning_of_time();
//        }
//    }
//
//    delete data;
//}
//
//TEST(perf_testing_openmp,max)
//{
//    perfOpenMP_general_test(&perf_testing_openmp_max);
//}
//
//double asum_testing_serial(const double data[], const int data_size)
//{
//    double asum = 0.0;
//    for (int i=0;i<data_size;i++)
//    {
//        asum+=std::abs(data[i]);
//    }
//    return asum;
//}
//
//double asum_testing_asumPrivates(const double data[], const int data_size)
//{
//    double shar_asum=0.0;
//    double priv_asum;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_asum)
//    {
//#endif
//        priv_asum=0.0;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (int i=0;i<data_size;i++)
//        {
//            priv_asum+=std::abs(data[i]);
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//        shar_asum+=priv_asum;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//    }
//#endif
//    return shar_asum;
//}
//
//double asum_testing_reduction(const double data[], const int data_size)
//{
//    double shar_asum=0.0;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) reduction(+:shar_asum)
//#endif
//    for (int i=0;i<data_size;i++)
//    {
//        shar_asum+=std::abs(data[i]);
//    }
//    return shar_asum;
//}
//
//double asum_testing_serial_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double asum=0.0;
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        asum+=stk::mesh::FortranBLAS<double>::asum((int)bucket_size,&data[i*bucket_size]);
//    }
//    return asum;
//}
//
//double asum_testing_asumPrivates_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_asum=0.0;
//    double priv_asum;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_asum)
//    {
//#endif
//        priv_asum=0.0;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (unsigned int i=0;i<num_buckets;i++)
//        {
//            priv_asum+=stk::mesh::FortranBLAS<double>::asum((int)bucket_size,&data[i*bucket_size]);
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//        shar_asum+=priv_asum;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//    }
//#endif
//    return shar_asum;
//}
//
//double asum_testing_reduction_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    double shar_asum=0.0;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) reduction(+:shar_asum)
//#endif
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        shar_asum+=stk::mesh::FortranBLAS<double>::asum((int)bucket_size,&data[i*bucket_size]);
//    }
//    return shar_asum;
//}
//
//void perf_testing_openmp_asum(const unsigned int num_buckets)
//{
//    std::cout<<std::endl<<"num_buckets,"<<num_buckets<<std::endl;
//    const double low_value = -2.07;
//    const double high_value = -4.27;
//    const unsigned int bucket_size = 512;
//    const int data_size = bucket_size*num_buckets;
//
//    const int num_numThreads = 4;
//    const int numThread_list [num_numThreads] = {1,10,20,30};
//
//    const unsigned int test_repeat = 750;
//
//    double * data = new double [data_size];
//    for (unsigned int i=0;i<(unsigned int)data_size;i++)
//    {
//        data[i]=low_value;
//    }
//
//    double expect_result = std::abs(low_value)*double(data_size);
//
//    const unsigned int seperation_increment = 100;
//
//    for (unsigned int i=seperation_increment;i<(unsigned int)(data_size-data_size%(int)seperation_increment);i+=seperation_increment)
//    {
//        data[i]=low_value+(high_value-low_value)*double(int(i)-(int)seperation_increment)/
//                double(-2*(int)seperation_increment+data_size-data_size%seperation_increment);
//        expect_result+=std::abs(data[i])-std::abs(low_value);
//    }
//
//    const double num_tol = expect_result*5.0e-9;
//    EXPECT_LT(std::abs(num_tol),std::abs(low_value*0.75));
//
//    unsigned int num_test_types = 3u;
//    TestFixture<double,double> test_fixture[num_test_types];
//    test_fixture[0].set_func_0_0(&asum_testing_serial);
//    test_fixture[0].set_func_0_1(&asum_testing_serial_buckets);
//    test_fixture[0].func_name=std::string ("serial");
//
//    test_fixture[1].set_func_0_0(&asum_testing_asumPrivates);
//    test_fixture[1].set_func_0_1(&asum_testing_asumPrivates_buckets);
//    test_fixture[1].func_name=std::string ("asumPrivates");
//
//    test_fixture[2].set_func_0_0(&asum_testing_reduction);
//    test_fixture[2].set_func_0_1(&asum_testing_reduction_buckets);
//    test_fixture[2].func_name=std::string ("reduction");
//
//    char coutBuffer [100];
//    sprintf(coutBuffer,"%15s , %10s , %10s , %10s , %10s , %10s\n","algorithm name","numThreads","NonBucket","NonBucket\%","Buckets","Buckets\%");
//    std::cout<<coutBuffer;
//
//    for (unsigned int j=0;j<num_test_types;j++)
//    {
//        double single_threaded_time_0 = 50.0;
//        double single_threaded_time_1 = 50.0;
//
//        for (int t=0;t<num_numThreads;t++)
//        {
//            if (j==0 && t>0) continue;
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            omp_set_num_threads(numThread_list[t]);
//#endif
//            for (unsigned int tr=0;tr<test_repeat;tr++)
//            {
//                EXPECT_NEAR(test_fixture[j].call_func_0(data,data_size),expect_result,num_tol);
//                EXPECT_NEAR(test_fixture[j].call_func_1(data,num_buckets,bucket_size),expect_result,num_tol);
//            }
//
//            if (t==0)
//            {
//                single_threaded_time_0=test_fixture[j].get_time(0);
//                single_threaded_time_1=test_fixture[j].get_time(1u);
//            }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            sprintf(coutBuffer,"%15s , %10d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#else
//            sprintf(coutBuffer,"%15s , %3d of %3d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    0,numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#endif
//            std::cout<<coutBuffer;
//
//            test_fixture[j].return_to_the_beginning_of_time();
//        }
//    }
//
//    delete[] data;
//}
//
//TEST(perf_testing_openmp,asum)
//{
//    perfOpenMP_general_test(&perf_testing_openmp_asum);
//}
//
//unsigned int imin_testing_serial(const double data[], const int data_size)
//{
//    unsigned int imin=0u;
//    double min=std::abs(data[0]);
//    double tmp;
//    for (int i=1;i<data_size;i++)
//    {
//        tmp=std::abs(data[i]);
//        if (tmp<min) {
//            imin=i;
//            min=tmp;
//        }
//    }
//    return imin;
//}
//
//unsigned int imin_testing_hiddenCritical(const double data[], const int data_size)
//{
//    unsigned int shar_imin=0u;
//    double shar_min=std::abs(data[0]);
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) private(priv_tmp)
//#endif
//    for (int i=1;i<data_size;i++)
//    {
//        priv_tmp=std::abs(data[i]);
//        if (priv_tmp<shar_min)
//        {
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//            if (priv_tmp<shar_min)
//            {
//#endif
//                shar_imin=i;
//                shar_min=priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            }
//#endif
//        }
//    }
//    return shar_imin;
//}
//
//unsigned int imin_testing_comparePrivates(const double data[], const int data_size)
//{
//    unsigned int shar_imin=0u;
//    double shar_min=std::abs(data[0]);
//    unsigned int priv_imin;
//    double priv_min;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_imin,priv_min,priv_tmp)
//    {
//#endif
//        priv_imin=0;
//        priv_min=std::abs(data[0]);
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (int i=1;i<data_size;i++)
//        {
//            priv_tmp=std::abs(data[i]);
//            if (priv_tmp<priv_min) {
//                priv_imin=i;
//                priv_min=priv_tmp;
//            }
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (priv_min<shar_min)
//        {
//#pragma omp critical
//            if (priv_min<shar_min)
//            {
//#endif
//                shar_imin=priv_imin;
//                shar_min=priv_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            }
//        }
//    }
//#endif
//    return shar_imin;
//}
//
//unsigned int imin_testing_reduction(const double data[], const int data_size)
//{
//    unsigned int shar_imin=0u;
//    double shar_min=std::abs(data[0]);
//    unsigned int priv_imin;
//    double priv_min;
//    double priv_tmp;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_tmp,priv_imin,priv_min)
//    {
//#endif
//        priv_imin=0;
//        priv_min=std::abs(data[0]);
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) reduction(min:shar_min)
//#endif
//        for (int i=0;i<data_size;i++)
//        {
//            priv_tmp=std::abs(data[i]);
//            if (priv_tmp<priv_min) {
//                priv_min=priv_tmp;
//                shar_min=priv_min;
//                priv_imin=i;
//            }
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (shar_min==priv_min)
//        {
//#endif
//            shar_imin=priv_imin;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        }
//    }
//#endif
//    return shar_imin;
//}
//
//unsigned int imin_testing_serial_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    unsigned int imin=0;
//    double min=std::abs(data[0]);
//    unsigned int tmp_imin;
//    double tmp_min;
//    for (unsigned int i=0u;i<num_buckets;i++)
//    {
//        tmp_imin=i*bucket_size+(unsigned int)stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size]);
//        tmp_min=std::abs(data[tmp_imin]);
//        if (tmp_min<min) {
//            min=tmp_min;
//            imin=tmp_imin;
//        }
//    }
//    return imin;
//}
//
//unsigned int imin_testing_hiddenCritical_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    unsigned int shar_imin=0;
//    double shar_min=std::abs(data[0]);
//    unsigned int priv_imin;
//    double priv_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel for schedule(static) private(priv_min,priv_imin)
//#endif
//    for (unsigned int i=0;i<num_buckets;i++)
//    {
//        priv_imin=i*bucket_size+(unsigned int)stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size]);
//        priv_min=std::abs(data[priv_imin]);
//        if (priv_min<shar_min)
//        {
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp critical
//#endif
//            if (priv_min<shar_min) {
//                shar_imin=priv_imin;
//                shar_min=priv_min;
//            }
//        }
//    }
//    return shar_imin;
//}
//
//unsigned int imin_testing_comparePrivates_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    unsigned int shar_imin=0;
//    double shar_min=std::abs(data[0]);
//    unsigned int priv_imin;
//    double priv_min;
//    unsigned int priv_tmp_imin;
//    double priv_tmp_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_imin,priv_min,priv_tmp_imin,priv_tmp_min)
//    {
//#endif
//        priv_imin=0;
//        priv_min=std::abs(data[0]);
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) nowait
//#endif
//        for (unsigned int i=0;i<num_buckets;i++)
//        {
//            priv_tmp_imin=i*bucket_size+(unsigned int)stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size]);
//            priv_tmp_min=std::abs(data[priv_tmp_imin]);
//            if (priv_tmp_min<priv_min)
//            {
//                priv_min=priv_tmp_min;
//                priv_imin=priv_tmp_imin;
//            }
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (priv_min<shar_min)
//        {
//#pragma omp critical
//            if (priv_min<shar_min)
//            {
//#endif
//                shar_min=priv_min;
//                shar_imin=priv_imin;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            }
//        }
//    }
//#endif
//    return shar_imin;
//}
//
//unsigned int imin_testing_reduction_buckets(const double data[], const unsigned int num_buckets, const unsigned int bucket_size)
//{
//    unsigned int shar_imin=0;
//    double shar_min=std::abs(data[0]);
//    unsigned int priv_imin;
//    double priv_min;
//    unsigned int priv_tmp_imin;
//    double priv_tmp_min;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp parallel private(priv_imin,priv_min,priv_tmp_imin,priv_tmp_min)
//    {
//#endif
//        priv_imin=0;
//        priv_min=std::abs(data[0]);
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//#pragma omp for schedule(static) reduction(min:shar_min)
//#endif
//        for (unsigned int i=0;i<num_buckets;i++)
//        {
//            priv_tmp_imin=i*bucket_size+(unsigned int)stk::mesh::FortranBLAS<double>::iamin((int)bucket_size,&data[i*bucket_size]);
//            priv_tmp_min=std::abs(data[priv_tmp_imin]);
//            if (priv_tmp_min<priv_min)
//            {
//                priv_min=priv_tmp_min;
//                priv_imin=priv_tmp_imin;
//                shar_min=priv_min;
//            }
//        }
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        if (priv_min==shar_min)
//        {
//#endif
//            shar_imin=priv_imin;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//        }
//    }
//#endif
//    return shar_imin;
//}
//
//void perf_testing_openmp_imin(const unsigned int num_buckets)
//{
//    std::cout<<std::endl<<"num_buckets,"<<num_buckets<<std::endl;
//    const double high_value = -4.27;
//    const double low_value = -3.14;
//    const unsigned int bucket_size = 512;
//    const int data_size = bucket_size*num_buckets;
//
//    const int num_numThreads = 4;
//    const int numThread_list [num_numThreads] = {1,10,20,30};
//
//    const unsigned int test_repeat = 750;
//
//    double * data = new double [data_size];
//    for (unsigned int i=0;i<(unsigned int)data_size;i++)
//    {
//        data[i]=high_value;
//    }
//
//    const unsigned int seperation_increment = 100;
//
//    for (unsigned int i=seperation_increment;i<(unsigned int)(data_size-data_size%(int)seperation_increment);i+=seperation_increment)
//    {
//        data[i]=high_value+(low_value-high_value)*double(int(i)-(int)seperation_increment)/
//                double(-2*(int)seperation_increment+data_size-data_size%seperation_increment);
//    }
//
//    unsigned int num_test_types = 4u;
//    TestFixture<double,unsigned int> test_fixture[num_test_types];
//    test_fixture[0].set_func_0_0(&imin_testing_serial);
//    test_fixture[0].set_func_0_1(&imin_testing_serial_buckets);
//    test_fixture[0].func_name=std::string ("serial");
//
//    test_fixture[1].set_func_0_0(&imin_testing_hiddenCritical);
//    test_fixture[1].set_func_0_1(&imin_testing_hiddenCritical_buckets);
//    test_fixture[1].func_name=std::string ("hiddenCritical");
//
//    test_fixture[2].set_func_0_0(&imin_testing_comparePrivates);
//    test_fixture[2].set_func_0_1(&imin_testing_comparePrivates_buckets);
//    test_fixture[2].func_name=std::string ("comparePrivates");
//
//    test_fixture[3].set_func_0_0(&imin_testing_reduction);
//    test_fixture[3].set_func_0_1(&imin_testing_reduction_buckets);
//    test_fixture[3].func_name=std::string ("reduction");
//
//    char coutBuffer [100];
//    sprintf(coutBuffer,"%15s , %10s , %10s , %10s , %10s , %10s\n","algorithm name","numThreads","NonBucket","NonBucket\%","Buckets","Buckets\%");
//    std::cout<<coutBuffer;
//
//    for (unsigned int j=0;j<num_test_types;j++)
//    {
//        double single_threaded_time_0 = 50.0;
//        double single_threaded_time_1 = 50.0;
//
//        for (int t=0;t<num_numThreads;t++)
//        {
//            if (j==0 && t>0) continue;
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            omp_set_num_threads(numThread_list[t]);
//#endif
//            for (unsigned int tr=0;tr<test_repeat;tr++)
//            {
//                EXPECT_EQ(data[test_fixture[j].call_func_0(data,data_size)],low_value);
//                EXPECT_EQ(data[test_fixture[j].call_func_1(data,num_buckets,bucket_size)],low_value);
//            }
//
//            if (t==0)
//            {
//                single_threaded_time_0=test_fixture[j].get_time(0);
//                single_threaded_time_1=test_fixture[j].get_time(1u);
//            }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMP_CPP
//            sprintf(coutBuffer,"%15s , %10d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#else
//            sprintf(coutBuffer,"%15s , %3d of %3d , %10f , %10f , %10f , %10f\n",
//                    test_fixture[j].func_name.data(),
//                    0,numThread_list[t],
//                    test_fixture[j].get_time(0),
//                    test_fixture[j].get_time(0)/single_threaded_time_0,
//                    test_fixture[j].get_time(1u),
//                    test_fixture[j].get_time(1u)/single_threaded_time_1);
//#endif
//            std::cout<<coutBuffer;
//
//            test_fixture[j].return_to_the_beginning_of_time();
//        }
//    }
//
//    delete[] data;
//}
//
//double omp_testing_asum(const int kappa,const double * x)
//{
//    double result=0.0;
//    for(int k=0;k<kappa;k++)
//    {
//        result+=std::abs(x[k]);
//    }
//    return result;
//}
//
//void omp_testing_axpy(const int kappa,const double alpha,const double * x,double * y)
//{
//    for(int k=0;k<kappa;k++)
//    {
//        y[k]+=alpha*x[k];
//    }
//}
//
//TEST(OMP,testing)
//{
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    std::cout<<"openmp verion: "<<_OPENMP<<std::endl;
//#endif
//
//    omp_set_dynamic(1);
//    omp_set_num_threads(omp_get_num_procs());
//#pragma omp parallel
//    {
//        if (omp_get_thread_num()==0) std::cout<<omp_get_num_threads()<<std::endl;
//    }
//    //    const int num_buckets = 100000;
//    //    const int buckets_size = 512;
//    //    const int data_size = num_buckets*buckets_size;
//    //
//    //    const int num_rep = 200;
//    //
//    //    const double d0 = 1.2;
//    //    const double alpha = 3.1;
//    //    const double d1 = 4.1;
//    //
//    //    double * data0 = new double [data_size];
//    //    double * data1 = new double [data_size];
//    //    for (int i=0;i<data_size;i++)
//    //    {
//    //        data0[i]=d0;
//    //        data1[i]=d1;
//    //    }
//    //
//    //    unsigned int num_hit_in_loop;
//    //
//    //#pragma omp parallel private(num_hit_in_loop)
//    //    {
//    //        num_hit_in_loop=0u;
//    //
//    //        for (int j=0;j<num_rep;j++)
//    //        {
//    //#pragma omp for schedule(static)
//    //            for (int i=0;i<num_buckets;i++)
//    //            {
//    //                omp_testing_axpy(buckets_size,alpha,&data0[i*buckets_size],&data1[i*buckets_size]);
//    //                num_hit_in_loop++;
//    //            }
//    //        }
//    //
//    //#pragma omp single
//    //        std::cout<<std::endl<<"axpy"<<std::endl;
//    //
//    //#pragma omp for schedule(static,1) ordered
//    //        for (int i=0;i<omp_get_num_threads();i++)
//    //        {
//    //#pragma omp ordered
//    //            {
//    //                std::cout<<omp_get_thread_num()<<" of "<<omp_get_num_threads()<<" ran "<<num_hit_in_loop<<std::endl;
//    //            }
//    //        }
//    //
//    //        if (omp_get_thread_num()==0)
//    //        {
//    //            std::cout<<"omp_get_num_threads "<<omp_get_num_threads()<<std::endl;
//    //            std::cout<<"omp_get_max_threads "<<omp_get_max_threads()<<std::endl;
//    //            std::cout<<"omp_get_num_procs "  <<omp_get_num_procs()  <<std::endl;
//    //        }
//    //    }
//    //
//    //#pragma omp parallel private(num_hit_in_loop)
//    //    {
//    //        num_hit_in_loop=0u;
//    //#pragma omp for schedule(static)
//    //        for (int i=0;i<num_buckets;i++)
//    //        {
//    //            for(int j=0;j<buckets_size;j++)
//    //            {
//    //                EXPECT_EQ(data0[i*buckets_size+j],d0);
//    //                EXPECT_NEAR(data1[i*buckets_size+j],d0*num_rep*alpha+d1,0.5);
//    //            }
//    //            num_hit_in_loop++;
//    //        }
//    //
//    //#pragma omp single
//    //        std::cout<<std::endl<<"validating"<<std::endl;
//    //
//    //#pragma omp for schedule(static,1) ordered
//    //        for (int i=0;i<omp_get_num_threads();i++)
//    //        {
//    //#pragma omp ordered
//    //            {
//    //                std::cout<<omp_get_thread_num()<<" of "<<omp_get_num_threads()<<" ran "<<num_hit_in_loop<<std::endl;
//    //            }
//    //        }
//    //
//    //        if (omp_get_thread_num()==0)
//    //        {
//    //            std::cout<<"omp_get_num_threads "<<omp_get_num_threads()<<std::endl;
//    //            std::cout<<"omp_get_max_threads "<<omp_get_max_threads()<<std::endl;
//    //            std::cout<<"omp_get_num_procs "  <<omp_get_num_procs()  <<std::endl;
//    //        }
//    //    }
//    //
//    //    double result=0.0;
//    //
//    //#pragma omp parallel private(num_hit_in_loop)
//    //    {
//    //        num_hit_in_loop=0u;
//    //#pragma omp for schedule(static) reduction(+:result)
//    //        for (int i=0;i<num_buckets;i++)
//    //        {
//    //            result+=omp_testing_asum(buckets_size,&data0[i*buckets_size]);
//    //            num_hit_in_loop++;
//    //        }
//    //
//    //#pragma omp single
//    //        std::cout<<std::endl<<"asum"<<std::endl;
//    //
//    //        EXPECT_NEAR(result,data_size*d0,0.5);
//    //
//    //#pragma omp for schedule(static,1) ordered
//    //        for (int i=0;i<omp_get_num_threads();i++)
//    //        {
//    //#pragma omp ordered
//    //            {
//    //                std::cout<<omp_get_thread_num()<<" of "<<omp_get_num_threads()<<" ran "<<num_hit_in_loop<<std::endl;
//    //            }
//    //        }
//    //
//    //        if (omp_get_thread_num()==0)
//    //        {
//    //            std::cout<<"omp_get_num_threads "<<omp_get_num_threads()<<std::endl;
//    //            std::cout<<"omp_get_max_threads "<<omp_get_max_threads()<<std::endl;
//    //            std::cout<<"omp_get_num_procs "  <<omp_get_num_procs()  <<std::endl;
//    //        }
//    //    }
//    //
//    //    delete[] data0;
//    //    delete[] data1;
//}
