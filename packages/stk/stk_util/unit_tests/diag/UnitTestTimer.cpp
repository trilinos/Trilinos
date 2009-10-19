#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <unistd.h>

#include <mpi.h>

// #include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/diag/Writer.hpp>

// using namespace use_case;


class UnitTestTimer : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE(UnitTestTimer);
  CPPUNIT_TEST(testUnit);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}
  
  void tearDown()
  {}

  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestTimer);

enum {
  TIMER_DOMAIN		= 0x00001000,		///< Enable domain timers
  TIMER_REGION		= 0x00002000,		///< Enable region timers
  TIMER_PROCEDURE	= 0x00004000,		///< Enable procedure timers
  TIMER_MECHANICS	= 0x00008000,		///< Enable mechanics timers
  TIMER_ALGORITHM	= 0x00010000,		///< Enable algorithm timers
  TIMER_SOLVER		= 0x00020000,		///< Enable solver timers
  TIMER_CONTACT		= 0x00040000,		///< Enable contact timers
  TIMER_MATERIAL	= 0x00080000,		///< Enable material timers
  TIMER_SEARCH		= 0x00100000,		///< Enable search timers
  TIMER_TRANSFER	= 0x00200000,		///< Enable transfer timers
  TIMER_ADAPTIVITY	= 0x00400000 		///< Enable adaptivity
};

enum {
  TIMER_UNUSED_1	= 0x00001000,		///< Enable unused 1
  TIMER_PROFILE_1	= 0x00002000,		///< Enable profile 1 timers
  TIMER_PROFILE_2	= 0x00004000,		///< Enable profile 2 timers
  TIMER_PROFILE_3	= 0x00008000,		///< Enable profile 3 timers
  TIMER_PROFILE_4	= 0x00010000,		///< Enable profile 4 timers
  TIMER_APP_1		= 0x00020000,		///< Enable application defined 1
  TIMER_APP_2		= 0x00040000,		///< Enable application defined 2
  TIMER_APP_3		= 0x00080000		///< Enable application defined 3
};

namespace {

double
quick_work()
{
  double x = 1.0;

  for (int i = 0; i < 10000; ++i) 
    x += std::sin((double) i);

  return x;
}


double
work()
{
  double x = 1.0;

  for (int i = 0; i < 100000; ++i) 
//  for (int i = 0; i < 100; ++i) 
    x += std::sin((double) i);

  return x;
}


stk::diag::TimerSet &
unitTestTimerSet()
{
  static stk::diag::TimerSet s_unitTestTimerSet(TIMER_REGION);

  return s_unitTestTimerSet;
}


stk::diag::TimerSet &
unitTestSecondTimerSet()
{
  static stk::diag::TimerSet s_unitTestSecondTimerSet(TIMER_APP_3);

  return s_unitTestSecondTimerSet;
}


stk::diag::Timer &unitTestTimer() {
  const std::string name("Unit test timer");
  static stk::diag::Timer s_unitTestTimer (stk::diag::createRootTimer(name, unitTestTimerSet()));

  return s_unitTestTimer;
}


struct RootObject
{
  RootObject()
    : m_timer("Root object", TIMER_REGION, unitTestTimer())
  {}

  stk::diag::Timer      m_timer;
};


struct Object 
{
  Object(const std::string &name, RootObject &root_object) 
    : m_id(0),
      m_name(name),
      m_timer(name, root_object.m_timer)
  {}
  
  Object(int id, const Object &parent) 
    : m_id(id),
      m_name(id_name(id)),
      m_timer(m_name, parent.m_timer)
  {}
  
  static std::string id_name(int id) {
    std::ostringstream s;
    s << "Object id " << id << " run";
    return s.str();
  }

  void run() {
    stk::diag::TimeBlock _time(m_timer);
    m_x += work();
  }
  
  int                   m_id;
  std::string           m_name;
  stk::diag::Timer      m_timer;
  double                m_x;
};

} // namespace <empty>


void
UnitTestTimer::testUnit()
{
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  stk::diag::TimeBlock root_time_block(unitTestTimer());

  std::ostringstream strout;
  
  // Create subtimer and test lap time
  {
    static stk::diag::Timer lap_timer("One second Wall time twice", unitTestTimer());
    
    stk::diag::TimeBlock _time(lap_timer);
    double x = quick_work();
    x = x;
    std::ostringstream oss;
    oss << x << std::endl;
    
    ::sleep(1);

    lap_timer.lap();
    
    stk::diag::MetricTraits<stk::diag::WallTime>::Type lap_time = lap_timer.getMetric<stk::diag::WallTime>().getLap();
  
    CPPUNIT_ASSERT(lap_time >= 1.0);

    ::sleep(1);

    lap_timer.stop();
    
    lap_time = lap_timer.getMetric<stk::diag::WallTime>().getLap();
  
    CPPUNIT_ASSERT(lap_time >= 2.0);
  }

  // 
  {
    static stk::diag::Timer run_timer("Run 100 times twice", unitTestTimer());
    
    for (int i = 0; i < 100; ++i) {
      stk::diag::TimeBlock _time(run_timer);
      work();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = run_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    CPPUNIT_ASSERT(lap_count == 100);
  }

  // Create second timer set
  {
    static stk::diag::Timer second_timer("Second timer set", unitTestTimer(), unitTestSecondTimerSet());
    static stk::diag::Timer second_timer_on_default("On default", second_timer);
    static stk::diag::Timer second_timer_on("On", TIMER_APP_3, second_timer);
    static stk::diag::Timer second_timer_off("Off", TIMER_APP_1, second_timer);
    
    stk::diag::TimeBlock _time(second_timer);
    stk::diag::TimeBlock _time1(second_timer_on_default);
    stk::diag::TimeBlock _time2(second_timer_on);
    stk::diag::TimeBlock _time3(second_timer_off);

    ::sleep(1);
  }

  // Grab previous subtimer and run 100 laps
  {
    static stk::diag::Timer run_timer("Run 100 times twice", unitTestTimer());
    
    for (int i = 0; i < 100; ++i) {
      stk::diag::TimeBlock _time(run_timer);
      work();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = run_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    CPPUNIT_ASSERT(lap_count == 200);
  }

  // Create root object
  RootObject root_object;
    
  {
    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = root_object.m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    CPPUNIT_ASSERT(lap_count == 0);
  }

  // Create object
  {
    Object time_object("One object", root_object);
    
    for (int i = 0; i < 100; ++i) {
      time_object.run();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = time_object.m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    CPPUNIT_ASSERT(lap_count == 100);
  }

  // Create object tree
  {
    std::vector<Object> object_vector;
    object_vector.push_back(Object("Object Tree", root_object));

    int id = 0;
    for (size_t i = 0; i < 2; ++i) {
      size_t ix = object_vector.size();
      object_vector.push_back(Object(id++, object_vector[0]));
      for (size_t j = 0; j < 2; ++j) {
        size_t jx = object_vector.size();
        object_vector.push_back(Object(id++, object_vector[ix]));
        for (int k = 0; k < 2; ++k) {    
          object_vector.push_back(Object(id++, object_vector[jx]));
        }
      }
    }
    
    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, false);
    
    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j) 
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    CPPUNIT_ASSERT_EQUAL(lap_count, stk::diag::MetricTraits<stk::diag::LapCount>::Type(0));

    for (size_t j = 0; j < object_vector.size(); ++j) 
      object_vector[j].run();

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, false);    

    lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j) 
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    CPPUNIT_ASSERT_EQUAL(lap_count, stk::diag::MetricTraits<stk::diag::LapCount>::Type(object_vector.size()));

    for (size_t i = 1; i < 100; ++i) 
      for (size_t j = 0; j < object_vector.size(); ++j) 
        object_vector[j].run();

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, false);    

    lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j) 
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);
  
    CPPUNIT_ASSERT_EQUAL(lap_count, stk::diag::MetricTraits<stk::diag::LapCount>::Type(100*object_vector.size()));

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, true);

    for (size_t i = 1; i < 100; ++i) 
      for (size_t j = 0; j < object_vector.size(); ++j) 
        object_vector[j].run();

    stk::diag::printTimersTable(strout, unitTestTimer(), stk::diag::METRICS_ALL, true);

    std::cout << strout.str() << std::endl;
    
//    dw().m(LOG_TIMER) << strout.str() << stk::diag::dendl;
  }
}
