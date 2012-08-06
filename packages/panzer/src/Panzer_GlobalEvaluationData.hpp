#ifndef __Panzer_GlobalEvaluationData_hpp__
#define __Panzer_GlobalEvaluationData_hpp__

namespace panzer {

/** This class is used by panzer to manage 
  * the data that is not contained in a workset.
  * It is often accessed by the gather/scatter
  * evaluators, where it is looked up by a string
  * identifier in the preEvaluate method. This lookup
  * is handled by the <code>GlobalEvaluatorDataContainer</code>.
  */
class GlobalEvaluationData {
public:
   virtual ~GlobalEvaluationData() = 0;

   virtual void ghostToGlobal(int mem) = 0;
   virtual void globalToGhost(int mem) = 0;

   virtual void initializeData() = 0;
};

/** Class that overides the communication primitives
  * to do nothing. This is used by the <code>LinearObjContainer</code>.
  */
class GlobalEvaluationData_Default : public GlobalEvaluationData {
public:
   virtual void ghostToGlobal(int mem) {}
   virtual void globalToGhost(int mem) {}
   virtual void initializeData() {}
};

}

#endif
