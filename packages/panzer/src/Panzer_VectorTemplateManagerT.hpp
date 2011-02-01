template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename BuilderT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::buildAndPushBackObjects(const BuilderT & builder)
{
   Teuchos::RCP<ManagerType> manager = Teuchos::rcp(new ManagerType);
   manager->buildObjects(builder);

   managers_.push_back(manager);
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::getAsBase(std::vector<Teuchos::RCP<BaseT> > & baseObjects)
{
   baseObjects.clear();

   // spin over the managers vector and extract the base objects
   typename std::vector<Teuchos::RCP<ManagerType> >::iterator itr;
   for(itr=managers_.begin();itr!=managers_.end();++itr) {
      Teuchos::RCP<BaseT> basePtr = (*itr)->getAsBase<EvalT>();
      baseObjects.push_back(basePtr);
   }
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::getAsBase(std::vector<Teuchos::RCP<const BaseT> > & baseObjects) const
{
   baseObjects.clear();

   // spin over the managers vector and extract the base objects
   typename std::vector<Teuchos::RCP<ManagerType> >::const_iterator itr;
   for(itr=managers_.begin();itr!=managers_.end();++itr) {
      Teuchos::RCP<const BaseT> basePtr = (*itr)->getAsBase<EvalT>();
      baseObjects.push_back(basePtr);
   }
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::
getAsObject(std::vector<Teuchos::RCP<typename boost::mpl::apply<ObjectT,EvalT>::type > > & objects)
{
   typedef typename boost::mpl::apply<ObjectT,EvalT>::type EvalObjectType;
   objects.clear();

   // spin over the managers vector and extract the base objects
   typename std::vector<Teuchos::RCP<ManagerType> >::iterator itr;
   for(itr=managers_.begin();itr!=managers_.end();++itr) {
      objects.push_back((*itr)->getAsObject<EvalT>());
   }
}

template <typename SeqTypes,typename BaseT,typename ObjectT>
template <typename EvalT>
void panzer::VectorTemplateManager<SeqTypes,BaseT,ObjectT>::
getAsObject(std::vector<Teuchos::RCP<const typename boost::mpl::apply<ObjectT,EvalT>::type > > & objects) const
{
   typedef typename boost::mpl::apply<ObjectT,EvalT>::type EvalObjectType;
   objects.clear();

   // spin over the managers vector and extract the base objects
   typename std::vector<Teuchos::RCP<ManagerType> >::const_iterator itr;
   for(itr=managers_.begin();itr!=managers_.end();++itr) {
      objects.push_back((*itr)->getAsObject<EvalT>());
   }
}
