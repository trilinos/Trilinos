template <typename SeqTypes,typename BaseT,typename ObjectT>
panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::VectorTemplateManager()
{
   int sz = Sacado::mpl::size<SeqTypes>::value;
   base_objects_.resize(sz);
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename BuilderT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::buildAndPushBackObjects(const BuilderT & builder)
{
   typedef PHX::TemplateManager<SeqTypes,BaseT,ObjectT> ManagerType;

   // get TemplateManager to build all local objects
   Teuchos::RCP<ManagerType> manager = Teuchos::rcp(new ManagerType);
   manager->buildObjects(builder);

   // copy local objects into my vectors
   std::size_t index=0;
   typename ManagerType::iterator itr = manager->begin();
   for(index=0;itr!=manager->end();itr++,index++)
      base_objects_[index].push_back(itr.rcp());
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::getAsBase(std::vector<Teuchos::RCP<BaseT> > & baseObjects)
{
   baseObjects.clear();

   // do a simple vector copy
   int idx = Sacado::mpl::find<SeqTypes,EvalT>::value;
   baseObjects = base_objects_[idx];
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::getAsBase(std::vector<Teuchos::RCP<const BaseT> > & baseObjects) const
{
   baseObjects.clear();

   int idx = Sacado::mpl::find<SeqTypes,EvalT>::value;
   const std::vector<Teuchos::RCP<BaseT> > & evalBase = base_objects_[idx];

   // do a simple vector copy
   typename std::vector<Teuchos::RCP<BaseT> >::const_iterator baseItr;
   for(baseItr=evalBase.begin();baseItr!=evalBase.end();++baseItr)
      baseObjects.push_back(*baseItr);
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::
getAsObject(std::vector<Teuchos::RCP<typename boost::mpl::apply<ObjectT,EvalT>::type > > & objects)
{
   typedef typename boost::mpl::apply<ObjectT,EvalT>::type EvalObjectType;
   objects.clear();

   // grab correctly index stl vector
   int idx = Sacado::mpl::find<SeqTypes,EvalT>::value;
   std::vector<Teuchos::RCP<BaseT> > & evalBase = base_objects_[idx];

   // loop over evalBase vector and cast to correct type
   typename std::vector<Teuchos::RCP<BaseT> >::iterator baseItr;
   for(baseItr=evalBase.begin();baseItr!=evalBase.end();++baseItr)
      objects.push_back(Teuchos::rcp_dynamic_cast<EvalObjectType>(*baseItr));
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::
getAsObject(std::vector<Teuchos::RCP<const typename boost::mpl::apply<ObjectT,EvalT>::type > > & objects) const
{
   typedef typename boost::mpl::apply<ObjectT,EvalT>::type EvalObjectType;
   objects.clear();

   // grab correctly index stl vector
   int idx = Sacado::mpl::find<SeqTypes,EvalT>::value;
   const std::vector<Teuchos::RCP<BaseT> > & evalBase = base_objects_[idx];

   // loop over evalBase vector and cast to correct type
   typename std::vector<Teuchos::RCP<BaseT> >::const_iterator baseItr;
   for(baseItr=evalBase.begin();baseItr!=evalBase.end();++baseItr)
      objects.push_back(Teuchos::rcp_dynamic_cast<const EvalObjectType>(*baseItr));
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
typename panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::iterator panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::begin()
{
   return panzer::VectorTemplateIterator<SeqTypes,BaseT,ObjectT>(*this,base_objects_.begin());
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
typename panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::iterator panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::end()
{
   return panzer::VectorTemplateIterator<SeqTypes,BaseT,ObjectT>(*this,base_objects_.end());
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
typename panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::const_iterator panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::begin() const
{
   return panzer::ConstVectorTemplateIterator<SeqTypes,BaseT,ObjectT>(*this,base_objects_.begin());
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
typename panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::const_iterator panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::end() const
{
   return panzer::ConstVectorTemplateIterator<SeqTypes,BaseT,ObjectT>(*this,base_objects_.end());
}
