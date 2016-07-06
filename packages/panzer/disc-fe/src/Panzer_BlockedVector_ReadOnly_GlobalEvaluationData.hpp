#ifndef __Panzer_BlockedVector_ReadOnly_GlobalEvaluationData_hpp__
#define __Panzer_BlockedVector_ReadOnly_GlobalEvaluationData_hpp__

#include "Panzer_ReadOnlyVector_GlobalEvaluationData.hpp"

namespace Thyra {
  // forward declaration
  template <typename> class DefaultProductVectorSpace;
}

namespace panzer {

/** This class encapsulates the needs of a gather operation to do a halo
  * exchange for blocked vectors.
  */
class BlockedVector_ReadOnly_GlobalEvaluationData : public ReadOnlyVector_GlobalEvaluationData {
public:

  BlockedVector_ReadOnly_GlobalEvaluationData();

  BlockedVector_ReadOnly_GlobalEvaluationData(const BlockedVector_ReadOnly_GlobalEvaluationData & src);

  BlockedVector_ReadOnly_GlobalEvaluationData(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ghostedSpace,
                                              const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > uniqueSpace,
                                              const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> > & gedBlocks);

  //! Virtual destructor
  virtual ~BlockedVector_ReadOnly_GlobalEvaluationData() {}

  /** Initialize this object using the sub GED objects. Also you must specify 
    * ghosted and unique spaces. At completion isInitialized_ will be set to true.
    *  
    * \param[in] ghostedSpace Must be a default product vector space defining the ghosted
    *                         vector
    * \param[in] ghostedSpace Must be a default product vector space defining the unique 
    *                         vector (currently ignored, included for future changes)
    * \param[in] gedBlocks GED objects that handle each block of the vector.
    */
  void initialize(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & ghostedSpace,
                  const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & uniqueSpace,
                  const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> > & gedBlocks);

  //! Is this object initialized
  virtual bool isInitialized() const { return isInitialized_; }

  /** For this class, this method does the halo exchange for the
    * vector.
    */
  virtual void globalToGhost(int mem);

  //! Initialize internal data for communication (clear the ghosted vector)
  virtual void initializeData();

  //! Set the unique vector
  virtual void setUniqueVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & uniqueVector);

  //! Get the unique vector
  virtual Teuchos::RCP<const Thyra::VectorBase<double> > getUniqueVector() const;

  //! Get the ghosted vector
  virtual Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const;

  //! How many blocks are in this GED
  size_t getBlockCount() const
  { return gedBlocks_.size(); }

  //! Get GED block (non const version)
  Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> getGEDBlock(int i)
  { return gedBlocks_[i]; }

  //! Get GED block (const version)
  Teuchos::RCP<const ReadOnlyVector_GlobalEvaluationData> getGEDBlock(int i) const
  { return gedBlocks_[i]; }

  //! No Dirichlet adjustment required
  bool requiresDirichletAdjustment() const { return false; }

private:

  bool isInitialized_;
 
  std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> > gedBlocks_;

  Teuchos::RCP<const Thyra::VectorBase<double> > uniqueVector_;

  Teuchos::RCP<const Thyra::DefaultProductVectorSpace<double> > ghostedSpace_;
};

}

#endif
