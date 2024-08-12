// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/markers/Marker.hpp>

#include <percept/PerceptMesh.hpp>
#include <adapt/BoundingRegion.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_io/StkMeshIoBroker.hpp>


#include <algorithm>
#include <functional>

#define DEBUG_PRINT 1
#define PRINTLN(a) do { if(DEBUG_PRINT && debug_print)  if (!bulkData_.parallel_rank()) std::cout << "  " << #a << " = " << a << std::endl; } while (0)
#define PRINTLN2(a,b) do { if(DEBUG_PRINT && debug_print)  if (!bulkData_.parallel_rank()) std::cout << "  " << #a << " = " << a << " " << #b << " = " << b << std::endl; } while (0)
#define PRINTLN3(a,b,c) do { if(DEBUG_PRINT && debug_print)  if (!bulkData_.parallel_rank()) std::cout << "  " << #a << " = " << a <<  " " << #b << " = " << b <<  " " << #c << " = " << c << std::endl; } while (0)

namespace percept{

//==========================================================================
// Class Definition
//==========================================================================
// Marker - computes refineField from errorIndicator field
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
  Marker::Marker(stk::mesh::BulkData& bulkData,
  MarkerInfo & markerInfo, stk::mesh::Selector *globalSelector) : bulkData_(bulkData),
                                                                  markerInfo_(markerInfo), doAdapt_(false),
                                                                  m_globalSelector(globalSelector)
{
  stk::mesh::MetaData & meta_data = bulkData.mesh_meta_data();
  if (!globalSelector)
    m_globalSelector = new stk::mesh::Selector(meta_data.universal_part());
  stk::mesh::EntityRank part_ranks[] = {stk::topology::ELEMENT_RANK, meta_data.side_rank()};
  for (unsigned irank=0; irank < 2; irank++)
    {
      std::ostringstream inactive_part_name;
      inactive_part_name << "refine_inactive_elements_part_" << static_cast<unsigned int>(part_ranks[irank]);
      //stk::mesh::Part* child_elements_part = m_eMesh.get_non_const_part(active_part_name);
      stk::mesh::Part* parent_elements_part = meta_data.get_part(inactive_part_name.str());
      if (!parent_elements_part) {
        //throw std::runtime_error("error - no parent_elements_part can be found");
        doAdapt_ = false;
      }
      else {
        adapterSelector_[part_ranks[irank]] = !stk::mesh::Selector(*parent_elements_part);
        doAdapt_ = true;
      }
    }
  if (!markerInfo_.transitionElementField_)
    {
      std::string name = "transition_element";
      if (meta_data.spatial_dimension() == 3)
        name = "transition_element_3";
      markerInfo_.transitionElementField_ = meta_data.get_field<percept::TransitionElementType::value_type>(stk::topology::ELEMENT_RANK, name);
    }
}


//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Marker::~Marker()
{
  delete m_globalSelector;
}

void Marker::setSelector(stk::mesh::Selector *sel)
{
  if (m_globalSelector) delete m_globalSelector;
  m_globalSelector = new stk::mesh::Selector(*sel);
}

stk::mesh::Selector *Marker::getSelector() { return m_globalSelector; }

  struct CompareErrIndRefFieldVec : public std::function<bool(ErrIndInfoTuple, ErrIndInfoTuple)> {
  bool operator()( ErrIndInfoTuple a,  ErrIndInfoTuple b)
  {
    return std::get<0>(a)  < std::get<0>(b);
  }
};

void
Marker::mark()
{
  if (!markerInfo_.useMarker_) {
    return;
  }
  zero_refine_field();
  computeMaxErrorIndicator();

  // define some common selectors;
  stk::mesh::Selector s_union =
    bulkData_.mesh_meta_data().locally_owned_part() &
    stk::mesh::selectField(*markerInfo_.errorIndicator_) & *m_globalSelector;

  std::vector<ErrIndInfoTuple > errIndRefFieldVec;

  stk::mesh::BucketVector const& elem_buckets =
    get_active_buckets(stk::topology::ELEMENT_RANK, s_union, false);

  const int ninfo = 5;
  int info[ninfo] = {0,0,0,0,0};

  // find refine field
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length = b.size();
    double * ei = stk::mesh::field_data(*markerInfo_.errorIndicator_, b);
    percept::RefineFieldType::value_type * refField = stk::mesh::field_data(*markerInfo_.refineField_, b);
    percept::RefineFieldType::value_type * refFieldOrig = stk::mesh::field_data(*markerInfo_.refineFieldOrig_, b);

    percept::TransitionElementType::value_type * transitionElement = 0;
    if (markerInfo_.transitionElementField_)
      transitionElement = stk::mesh::field_data(*markerInfo_.transitionElementField_, b);

    int *refineLevelField = stk::mesh::field_data(*markerInfo_.refineLevelField_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      int *refFieldVec = 0;
      bool tooRefined = false;

      if (refine_element(ei[k])) {
        ++info[0];

          refFieldOrig[k] = 1;
          if (refineLevelField && refineLevelField[k] >= markerInfo_.maxRefinementLevel_) {
            ++info[1];
            if (transitionElement && !transitionElement[k])
              {
                ++info[2];
                refField[k] = 0;
                tooRefined = true;
              }
            else
              {
                ++info[3];
                refField[k] = 1;
                tooRefined = false;
              }
          }
          else {
            ++info[4];
            refField[k] = 1;
          }
        }
      else if (unrefine_element(ei[k])) {
          refFieldOrig[k] = -1;
          refField[k] = -1;
        }
      else {
        refFieldOrig[k] = 0;
        refField[k] = 0;
      }

      // limit by geometric box
      if (!(markerInfo_.boundingRegion_.is_null()) &&
          !(markerInfo_.boundingRegion_->withinGeometricLimits(b[k]))) {
        refFieldOrig[k] = 0;
        refField[k] = 0;
        tooRefined = true;
      }

      refFieldVec = &refField[k];

      errIndRefFieldVec.push_back(ErrIndInfoTuple(&ei[k], refFieldVec, tooRefined));
    }
  }

  checkElementBudget(errIndRefFieldVec);


  if (markerInfo_.debug_)
    {
      int g_info[ninfo];
      stk::all_reduce_sum(bulkData_.parallel(), &info[0], &g_info[0], ninfo);
      if (bulkData_.parallel_rank() == 0)
        {
          std::cout << "Marker:: info= refElem tooRef +isTE -isTE !tooRef: "
                << g_info[0] << ", "
                << g_info[1] << ", "
                << g_info[2] << ", "
                << g_info[3] << ", "
                << g_info[4]
                << std::endl;
        }
    }

  // ensure parallel consistency
  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(markerInfo_.refineField_);
    fields.push_back(markerInfo_.refineFieldOrig_);
    stk::mesh::communicate_field_data(bulkData_.aura_ghosting(), fields);
  }

}

void
Marker::computeMaxErrorIndicator()
{
  // ensure parallel consistency
  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(markerInfo_.errorIndicator_);
    stk::mesh::communicate_field_data(bulkData_.aura_ghosting(), fields);
  }

  // find max to use for normalization - don't actually normalize since this
  // makes it difficult to visualize time-sequences of the error indicator
  double maxErrorIndicator = -std::numeric_limits<double>::max();
  double minErrorIndicator = +std::numeric_limits<double>::max();

  stk::mesh::MetaData & meta_data = bulkData_.mesh_meta_data();

  // define some common selectors;
  stk::mesh::Selector s_union =
    meta_data.locally_owned_part()
    & stk::mesh::selectField(*markerInfo_.errorIndicator_) & *m_globalSelector;

  bool get_all_including_parents = false;
  stk::mesh::BucketVector const& prior_elem_buckets =
    get_active_buckets( stk::topology::ELEMENT_RANK, s_union, get_all_including_parents );

  for ( stk::mesh::BucketVector::const_iterator ib = prior_elem_buckets.begin() ;
        ib != prior_elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * ei = stk::mesh::field_data(*markerInfo_.errorIndicator_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      maxErrorIndicator = std::max(ei[k], maxErrorIndicator);
      minErrorIndicator = std::min(ei[k], minErrorIndicator);
    }
  }

  stk::all_reduce_max(bulkData_.parallel(), &maxErrorIndicator, &markerInfo_.maxErrorIndicator_, 1);
  stk::all_reduce_min(bulkData_.parallel(), &minErrorIndicator, &markerInfo_.minErrorIndicator_, 1);

  if (markerInfo_.debug_)
    {
      if (!bulkData_.parallel_rank()) std::cout << "minErrorIndicator= " << markerInfo_.minErrorIndicator_ << std::endl;
      if (!bulkData_.parallel_rank()) std::cout << "maxErrorIndicator= " << markerInfo_.maxErrorIndicator_ << std::endl;
    }
}


static void printNref(stk::mesh::BulkData& bulkData_, const std::string& msg)
{
  stk::mesh::MetaData & meta_data = bulkData_.mesh_meta_data();
  percept::PerceptMesh eMesh(&meta_data, &bulkData_, true);
  percept::RefineFieldType *refine_field = eMesh.get_fem_meta_data()->  get_field<percept::RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
  stk::mesh::Selector on_locally_owned_part =  meta_data.locally_owned_part();
  const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );

  size_t nref=0;
  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      if (on_locally_owned_part(**k))
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];

              if (eMesh.numChildren(element))
                continue;

              //int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
              int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
              if (refine_field_elem[0] >= 1)
                {
                  ++nref;
                }
            }
        }
    }
  stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();

  stk::all_reduce( pm, stk::ReduceSum<1>( &nref ) );
  if (eMesh.get_rank()==0)
    std::cout << "P[" << eMesh.get_rank() << "] printNref nref1 for this iter= " << nref << " msg= " << msg << std::endl;
}

struct IterationInfo
{
  int iter;
  double currentErrorThreshold;
  double currentErrorThresholdMin;
  double currentErrorThresholdMax;
  size_t numCurrentElements;
  size_t numElementLimit;
  int convCase;
  IterationInfo(int iter, double c0, double c1, double c2, size_t n0, size_t n1, int convCase) :
    iter(iter), currentErrorThreshold(c0), currentErrorThresholdMin(c1), currentErrorThresholdMax(c2),
    numCurrentElements(n0), numElementLimit(n1), convCase(convCase)
  {
  }

};

static int64_t my_abs(int64_t a)
{
  if (a < 0) return -a;
  return a;
}

// checks and enforces any specified element budget
void Marker::checkElementBudget(std::vector<ErrIndInfoTuple >& errIndRefFieldVec)
{
  std::vector<IterationInfo> info;
  stk::mesh::MetaData & meta_data = bulkData_.mesh_meta_data();
  percept::PerceptMesh eMesh(&meta_data, &bulkData_, true);


  // FIXME - if frac is == 1.0, we need to restrict to no refinement, not just return
  // check for refining too many elements
  const double frac = markerInfo_.maxRefinementNumberOfElementsFraction_ ;
  if (frac <= 1.0) return;

  const bool debug_print = markerInfo_.debug_;

  if (debug_print) printNref(bulkData_, "before");

  // approximate estimate of elements refined from an individual element
  const unsigned numRefinedElemPerElem = meta_data.spatial_dimension() == 3 ? 8 : 4;

  int64_t numCurrentElements = count_elements(meta_data.locally_owned_part(), false);

      // not strictly necessary, but, could be useful for performance improvement in future
      //   (allows early breaks from some loops, but sort itself may be more expensive)
      //CompareErrIndRefFieldVec comp;
      //std::sort(errIndRefFieldVec.begin(), errIndRefFieldVec.end(), comp);

      // binary search (note: using iteration avoids the need for a global sort, at the
      //   expense of multiple parallel reductions, but only of a single scalar integer)
      // double currentErrorThresholdMax = R;
      // double currentErrorThresholdMin = 0.0;
      // double currentErrorThreshold = currentErrorThresholdMax;

  double currentErrorThresholdMax = 0.0;
  double currentErrorThresholdMin = 0.0;
  getThresholdMinMax(currentErrorThresholdMin, currentErrorThresholdMax);
  double currentErrorThreshold = currentErrorThresholdMin;

  PRINTLN2(currentErrorThresholdMin,currentErrorThresholdMax);

  int64_t numRefinedGlobal = estimateNewElements(currentErrorThreshold, errIndRefFieldVec);
  int64_t numTotalNewElements = numRefinedGlobal*numRefinedElemPerElem + numCurrentElements;
  const int64_t numTotalElementLimit = (int64_t)((frac)*double(markerInfo_.numInitialElements_));

  info.push_back(IterationInfo(-1, currentErrorThreshold, currentErrorThresholdMin, currentErrorThresholdMax, numTotalNewElements, numTotalElementLimit, 0));

  if (markerInfo_.debug_)
    {
      int64_t numRefinedGlobalMax = estimateNewElements(currentErrorThresholdMax, errIndRefFieldVec);
      int64_t numTotalNewElementsMax = numRefinedGlobalMax*numRefinedElemPerElem + numCurrentElements;
      if (!bulkData_.parallel_rank())
        std::cout
          << "Adapt: Predicted number of new elements= " << numTotalNewElements << " (max=  " << numTotalNewElementsMax
          << " currentErrorThresholdMax= " << currentErrorThresholdMax
          <<  " numRefinedGlobalMax= " << numRefinedGlobalMax
          << " numCurrentElements= " << numCurrentElements
          << "\n        Target (budget) number of new elements = " << numTotalElementLimit << std::endl;
    }
  PRINTLN2(numTotalNewElements,numTotalElementLimit);

  // too many elements already - skip adjustments to meet element budget
  if (numCurrentElements >= numTotalElementLimit)
    {
      markUsing(currentErrorThreshold, errIndRefFieldVec, false);
      if (!bulkData_.parallel_rank())
        std::cout << "Adapt: WARNING: element budget reached before iteration, no refine will be done, numCurrentElements= " << numCurrentElements
                  << " numTotalElementLimit= " << numTotalElementLimit << std::endl;
    }
  // estimated new elements exceeds limit on total elements,
  // go ahead and try to adjust to the threshold to have fewer elements
  else if (numTotalNewElements > numTotalElementLimit)
    {
      unsigned numIter = 0;
      const unsigned maxIter = markerInfo_.maxMarkerIterations_;
      bool found = false;
      std::vector<int64_t> numTotalElementsHistory;
      while(true)
        {
          currentErrorThreshold = 0.5*(currentErrorThresholdMin + currentErrorThresholdMax);
          PRINTLN3(currentErrorThresholdMin, currentErrorThresholdMax, currentErrorThreshold);

          numRefinedGlobal = estimateNewElements(currentErrorThreshold, errIndRefFieldVec);
          numTotalNewElements = numRefinedGlobal*numRefinedElemPerElem + numCurrentElements;
          numTotalElementsHistory.push_back(numTotalNewElements);

          if (markerInfo_.debug_)
            {
              if (!bulkData_.parallel_rank()) std::cout << "Adapt: restricting refined elements... iter= " << numIter
                                                        << " currentErrorThreshold= " << currentErrorThreshold 
                                                        << " numCurrentElements= " << numCurrentElements
                                                        << " numRefinedGlobal= " << numRefinedGlobal
                                                        << " numTotalNewElements= " << numTotalNewElements
                                                        << " numTotalElementLimit= " << numTotalElementLimit
                                                        << " numInitialElements_=  " << markerInfo_.numInitialElements_
                                                        << std::endl;
            }


          PRINTLN3(numIter, numTotalNewElements, numTotalElementLimit);
          int caseInfo = 0;
          if (markerInfo_.debug_ && eMesh.get_rank() == 0)
            {
              std::cout << "iter= " << numIter << " numTotalNewElements= " << numTotalNewElements << " numTotalElementLimit= " << numTotalElementLimit
                        << "\n (my_abs(numTotalNewElements - numTotalElementLimit) = " << (my_abs(numTotalNewElements - numTotalElementLimit))
                        << "\n numTotalNewElements == numTotalElementLimit = " << (numTotalNewElements == numTotalElementLimit)
                        << "\n numTotalNewElements < (int64_t)(1.05*double(numTotalElementLimit)) = " << (numTotalNewElements < (int64_t)(1.05*double(numTotalElementLimit)))
                        << "\n my_abs(numTotalNewElements - numTotalElementLimit) < numRefinedElemPerElem)= " << (my_abs(numTotalNewElements - numTotalElementLimit) < numRefinedElemPerElem)
                        << std::endl;
            }

          if (   numTotalNewElements == numTotalElementLimit
              || (numIter >= maxIter && numTotalNewElements < (int64_t)(1.05*double(numTotalElementLimit)))
                 || my_abs(numTotalNewElements - numTotalElementLimit) < numRefinedElemPerElem)
            {
              if (debug_print)
                if (!bulkData_.parallel_rank()) std::cout << "Adapt: converged, breaking (numTotalNewElements == numTotalElementLimit) ... " << std::endl;
              found = true;
              caseInfo = 3;
            }
          // too many elements - raise the minimum threshold
          else if (numTotalNewElements > numTotalElementLimit)
            {
              currentErrorThresholdMin = currentErrorThreshold;
              caseInfo = 1;
            }
          // too few elements - lower the max threshold
          else if (numTotalNewElements < numTotalElementLimit)
            {
              currentErrorThresholdMax = currentErrorThreshold;
              caseInfo = 2;
            }
          VERIFY_OP_ON(caseInfo, !=, 0, "bad caseInfo");
          info.push_back(IterationInfo(numIter, currentErrorThreshold, currentErrorThresholdMin, currentErrorThresholdMax, numTotalNewElements, numTotalElementLimit, caseInfo));
          if (found)
            break;

          // check if we can exit the while loop
          if (numIter > maxIter)
            {
              if (!bulkData_.parallel_rank()) {
                std::cout << "Adapt: WARNING: too many iterations in marking, no refine will be done, numTotalNewElementsBest= "
                          << numTotalNewElements << " numTotalElementLimit= " << numTotalElementLimit << std::endl;
                std::cout << "  History of total number of elements: ";
                for (unsigned i=0; i<numTotalElementsHistory.size(); i++) {
                  std::cout << numTotalElementsHistory[i] << " ";
                }
                std::cout << std::endl;
              }
              break;
            }
          ++numIter;
        } // while(true)

      PRINTLN2(currentErrorThreshold, found);
      markUsing(currentErrorThreshold, errIndRefFieldVec, found);
    }

  if (eMesh.get_rank() == 0)
    {
      int s = 25;
      int p = 20;
      std::ostringstream out;
      out << "Marker binary search convergence history:\n"
          << std::setprecision(p)
          << std::right
          << std::setw(s) << "iter"
          << std::setw(s) << "current error"
          << std::setw(s) << "min error"
          << std::setw(s) << "max error"
          << std::setw(s) << "# elements"
          << std::setw(s) << "element limit"
          << std::setw(s) << "conv type"
          << std::endl;
      for (unsigned ii=0; ii < info.size(); ++ii)
        {
          out << std::right
              << std::setw(s) << info[ii].iter
              << std::setw(s) << info[ii].currentErrorThreshold
              << std::setw(s) << info[ii].currentErrorThresholdMin
              << std::setw(s) << info[ii].currentErrorThresholdMax
              << std::setw(s) << info[ii].numCurrentElements
              << std::setw(s) << info[ii].numElementLimit
              << std::setw(s) << info[ii].convCase
              << std::endl;
        }
      std::cout << out.str() << std::endl;
    }

  if (debug_print) printNref(bulkData_, "after");
}

stk::mesh::BucketVector const& Marker::get_active_buckets( stk::mesh::EntityRank rank,
                                                           const stk::mesh::Selector & selector ,
                                                           bool get_all) const
{
  stk::mesh::BulkData const& mesh = bulkData_;

  if (!get_all && doAdapt_)
    {
      stk::mesh::Selector new_selector = selector;
      if (rank != stk::topology::NODE_RANK)
        {
          // adapterSelector_ avoids parent elements
          new_selector = selector & adapterSelector_[rank];
        }
      return mesh.get_buckets(rank, new_selector);
    }
  else
    {
      return mesh.get_buckets(rank, selector);
    }
}

size_t Marker::count_elements(stk::mesh::Selector sel, bool count_all_including_parents )
{
  stk::mesh::BucketVector const& elem_buckets =
    //bulkData_.get_buckets(stk::topology::ELEMENT_RANK, s_union);
    get_active_buckets(stk::topology::ELEMENT_RANK, sel, count_all_including_parents);

  size_t count = 0;
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length = b.size();
    count += length;
  }
  size_t g_count = 0;
  stk::all_reduce_sum(bulkData_.parallel(), &count, &g_count, 1);

  return g_count;
}

void Marker::zero_refine_field()
{

  if (!markerInfo_.useMarker_)
    {
      return;
    }
  //stk::mesh::MetaData & meta_data = bulkData_.mesh_meta_data();

  // define some common selectors;
  stk::mesh::Selector s_union =
    //meta_data.locally_owned_part() &
    stk::mesh::selectField(*markerInfo_.refineField_) & *m_globalSelector;

  bool get_all = true;
  stk::mesh::BucketVector const& elem_buckets =
    //bulkData_.get_buckets(stk::topology::ELEMENT_RANK, s_union);
    get_active_buckets(stk::topology::ELEMENT_RANK, s_union, get_all);

  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin() ;
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length = b.size();
    percept::RefineFieldType::value_type * refField = 0;
    percept::RefineFieldType::value_type * refFieldOrig = 0;
    refField = stk::mesh::field_data(*markerInfo_.refineField_, b);
    refFieldOrig = stk::mesh::field_data(*markerInfo_.refineFieldOrig_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      refField[k] = 0;
      refFieldOrig[k] = 0;
    }
  }

}

} // namespace percept
