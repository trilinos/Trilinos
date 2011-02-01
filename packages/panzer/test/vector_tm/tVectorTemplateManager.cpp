#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_VectorTemplateManager.hpp"
#include "Sacado_mpl_vector.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

struct TypeA  { static std::string name() { return "ambitious"; } };
struct TypeB  { static std::string name() { return "relaxed"; } };
struct TypeAB { static std::string name() { return "confused"; } };

typedef Sacado::mpl::vector<TypeA, TypeB, TypeAB> EvalTypes;

class BaseClass { public: 
   virtual ~BaseClass() {} 
   virtual std::string getName() const = 0;
};

template <typename EvalT>
class SampleType : public BaseClass {
public:
   virtual std::string getName() const = 0;
};

template <typename EvalT>
class SampleTypeA : public SampleType<EvalT> {
   std::string extra_;
public:
   SampleTypeA(const std::string & s) : extra_(s) {}
   virtual std::string getName() const { return "SampleA::"+extra_+"_"+EvalT::name(); } 
};

template <typename EvalT>
class SampleTypeB : public SampleType<EvalT> {
   std::string extra_;
public:
   SampleTypeB(const std::string & s) : extra_(s) {}
   virtual std::string getName() const { return "SampleB::"+extra_+"_"+EvalT::name(); } 
};

struct SampleTypeA_Builder {
   std::string extra_;
   SampleTypeA_Builder(std::string s) : extra_(s) {}

   template <typename EvalT>
   Teuchos::RCP<BaseClass> build() const 
   { return Teuchos::rcp(new SampleTypeA<EvalT>(extra_)); }
};

struct SampleTypeB_Builder {
   std::string extra_;
   SampleTypeB_Builder(std::string s) : extra_(s) {}

   template <typename EvalT>
   Teuchos::RCP<BaseClass> build() const 
   { return Teuchos::rcp(new SampleTypeB<EvalT>(extra_)); }
};

void getTruthVector(std::string evalName,std::vector<std::string> & v)
{
  v.clear();

  v.push_back("SampleA::dog_"+evalName); 
  v.push_back("SampleA::cat_"+evalName); 
  v.push_back("SampleB::horse_"+evalName); 
  v.push_back("SampleB::mouse_"+evalName); 

  std::sort(v.begin(),v.end());
}

template <typename ClassT>
void getVectorUnderTest(const std::vector<ClassT> & in,std::vector<std::string> & out)
{
   out.clear();
   for(std::size_t i=0;i<in.size();i++)
      out.push_back(in[i]->getName());

   std::sort(out.begin(),out.end());
}

TEUCHOS_UNIT_TEST(tVectorTemplateManager, test)
{
   // This test is kind icky.  The point here is that each
   // SampleType*<EvalT> returns a unique string. We build all of these
   // objects and verify that their strings are correct. This basically
   // doesn't check the order things are returned in but sorts the vectors
   // so they can be easily compared.

   typedef panzer::VectorTemplateManager<EvalTypes,BaseClass,SampleType<_> > VectorTM;

   VectorTM vectorTM;
   const VectorTM & constVectorTM = vectorTM;

   vectorTM.buildAndPushBackObjects(SampleTypeA_Builder("dog"));
   vectorTM.buildAndPushBackObjects(SampleTypeA_Builder("cat"));
   vectorTM.buildAndPushBackObjects(SampleTypeB_Builder("horse"));
   vectorTM.buildAndPushBackObjects(SampleTypeB_Builder("mouse"));

   std::vector<std::string> truthVector;
   std::vector<std::string> testVector;
   {
      std::vector<Teuchos::RCP<BaseClass> > baseObjs;

      // test AB
      vectorTM.getAsBase<TypeAB>(baseObjs);
      getVectorUnderTest(baseObjs,testVector);
      getTruthVector(TypeAB::name(),truthVector);

      TEST_EQUALITY(truthVector.size(),testVector.size());
      TEST_ASSERT(std::equal(truthVector.begin(),truthVector.end(),
                             testVector.begin()));

      // test A
      vectorTM.getAsBase<TypeA>(baseObjs);
      getVectorUnderTest(baseObjs,testVector);
      getTruthVector(TypeA::name(),truthVector);

      TEST_EQUALITY(truthVector.size(),testVector.size());
      TEST_ASSERT(std::equal(truthVector.begin(),truthVector.end(),
                             testVector.begin()));
   }

   {
      std::vector<Teuchos::RCP<const BaseClass> > baseObjs;

      // test AB
      constVectorTM.getAsBase<TypeAB>(baseObjs);
      getVectorUnderTest(baseObjs,testVector);
      getTruthVector(TypeAB::name(),truthVector);

      TEST_EQUALITY(truthVector.size(),testVector.size());
      TEST_ASSERT(std::equal(truthVector.begin(),truthVector.end(),
                             testVector.begin()));

      // test A
      constVectorTM.getAsBase<TypeA>(baseObjs);
      getVectorUnderTest(baseObjs,testVector);
      getTruthVector(TypeA::name(),truthVector);

      TEST_EQUALITY(truthVector.size(),testVector.size());
      TEST_ASSERT(std::equal(truthVector.begin(),truthVector.end(),
                             testVector.begin()));
   }

   {
      std::vector<Teuchos::RCP<SampleType<TypeB> > > baseObjs;

      // test B
      vectorTM.getAsObject<TypeB>(baseObjs);
      getVectorUnderTest(baseObjs,testVector);
      getTruthVector(TypeB::name(),truthVector);

      TEST_EQUALITY(truthVector.size(),testVector.size());
      TEST_ASSERT(std::equal(truthVector.begin(),truthVector.end(),
                             testVector.begin()));
   }

   {
      std::vector<Teuchos::RCP<const SampleType<TypeA> > > baseObjs;

      // test A
      constVectorTM.getAsObject<TypeA>(baseObjs);
      getVectorUnderTest(baseObjs,testVector);
      getTruthVector(TypeA::name(),truthVector);

      TEST_EQUALITY(truthVector.size(),testVector.size());
      TEST_ASSERT(std::equal(truthVector.begin(),truthVector.end(),
                             testVector.begin()));
   }
}

}
