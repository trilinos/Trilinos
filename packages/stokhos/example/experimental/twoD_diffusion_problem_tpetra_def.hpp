// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
twoD_diffusion_problem(
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
  LocalOrdinal n, LocalOrdinal d, 
  BasisScalar s, BasisScalar mu, 
  bool log_normal_,
  bool eliminate_bcs_) :
  mesh(n*n),
  log_normal(log_normal_),
  eliminate_bcs(eliminate_bcs_)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::global_size_t;

  //////////////////////////////////////////////////////////////////////////////
  // Construct the mesh.  
  // The mesh is uniform and the nodes are numbered
  // LEFT to RIGHT, DOWN to UP.
  //
  // 5-6-7-8-9
  // | | | | |
  // 0-1-2-3-4
  /////////////////////////////////////////////////////////////////////////////
  MeshScalar xyLeft = -.5;
  MeshScalar xyRight = .5;
  h = (xyRight - xyLeft)/((MeshScalar)(n-1));
  Array<GlobalOrdinal> global_dof_indices;
  for (GlobalOrdinal j=0; j<n; j++) {
    MeshScalar y = xyLeft + j*h;
    for (GlobalOrdinal i=0; i<n; i++) {
      MeshScalar x = xyLeft + i*h;
      GlobalOrdinal idx = j*n+i;
      mesh[idx].x = x;
      mesh[idx].y = y;
      if (i == 0 || i == n-1 || j == 0 || j == n-1)
	mesh[idx].boundary = true;
      if (i != 0)
	mesh[idx].left = idx-1;
      if (i != n-1)
	mesh[idx].right = idx+1;
      if (j != 0)
	mesh[idx].down = idx-n;
      if (j != n-1)
	mesh[idx].up = idx+n;
      if (!(eliminate_bcs && mesh[idx].boundary))
	global_dof_indices.push_back(idx);
    }
  }
  
  // Solution vector map
  global_size_t n_global_dof = global_dof_indices.size();
  int n_proc = comm->getSize();
  int proc_id = comm->getRank();
  size_t n_my_dof = n_global_dof / n_proc;
  if (proc_id == n_proc-1)
    n_my_dof += n_global_dof % n_proc;
  ArrayView<GlobalOrdinal> my_dof = 
    global_dof_indices.view(proc_id*(n_global_dof / n_proc), n_my_dof);
  x_map = Tpetra::createNonContigMap<LocalOrdinal,GlobalOrdinal>(my_dof, comm);

  // Initial guess, initialized to 0.0
  x_init = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(x_map);
  x_init->putScalar(0.0);

  // Parameter vector map
  p_map = Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(d, comm);

  // Response vector map
  g_map = Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(1, comm);

  // Initial parameters
  p_init = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(p_map);
  p_init->putScalar(0.0);

  // Parameter names
  p_names = Teuchos::rcp(new Array<std::string>(d));
  for (LocalOrdinal i=0;i<d;i++) {
    std::stringstream ss;
    ss << "KL Random Variable " << i+1;
    (*p_names)[i] = ss.str(); 
  }

  // Build Jacobian graph
  size_t NumMyElements = x_map->getLocalNumElements();
  ArrayView<const GlobalOrdinal> MyGlobalElements = 
    x_map->getLocalElementList ();
  graph = rcp(new Tpetra_CrsGraph(x_map, 5));
  for (size_t i=0; i<NumMyElements; ++i ) {

    // Center
    GlobalOrdinal global_idx = MyGlobalElements[i];
    graph->insertGlobalIndices(global_idx, arrayView(&global_idx, 1));

    if (!mesh[global_idx].boundary) {
      // Down
      if (!(eliminate_bcs && mesh[mesh[global_idx].down].boundary))
	graph->insertGlobalIndices(global_idx, 
				   arrayView(&mesh[global_idx].down,1));

      // Left
      if (!(eliminate_bcs && mesh[mesh[global_idx].left].boundary))
	graph->insertGlobalIndices(global_idx, 
				   arrayView(&mesh[global_idx].left,1));

      // Right
      if (!(eliminate_bcs && mesh[mesh[global_idx].right].boundary))
	graph->insertGlobalIndices(global_idx, 
				   arrayView(&mesh[global_idx].right,1));

      // Up
      if (!(eliminate_bcs && mesh[mesh[global_idx].up].boundary))
	graph->insertGlobalIndices(global_idx, 
				   arrayView(&mesh[global_idx].up,1));
    }
  }
  graph->fillComplete();

  // Construct deterministic operator
  A = rcp(new Tpetra_CrsMatrix(graph));
 
  // Construct the RHS vector.
  b = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(x_map);
  ArrayRCP<Scalar> b_view = b->get1dViewNonConst();
  for(size_t i=0; i<NumMyElements; ++i) {
    GlobalOrdinal global_idx = MyGlobalElements[i];
    if (mesh[global_idx].boundary)
      b_view[i] = 0;
    else 
      b_view[i] = 1;
  }

  // Diffusion functions
  klFunc = rcp(new KL_Diffusion_Func(xyLeft, xyRight, mu, s, 1.0, d));
  lnFunc = rcp(new LogNormal_Diffusion_Func<KL_Diffusion_Func>(*klFunc));
}

// Overridden from EpetraExt::ModelEvaluator
template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const 
	     typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_Map>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_x_map() const
{
  return x_map;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const 
	     typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_Map>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_f_map() const
{
  return x_map;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const 
	     typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_Map>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_p_map(LocalOrdinal l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_problem::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_map;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const 
	     typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_Map>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_g_map(LocalOrdinal l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_problem::get_g_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return g_map;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Teuchos::Array<std::string> >
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_p_names(LocalOrdinal l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_problem::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return p_names;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const 
	     typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_Vector>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_x_init() const
{
  return x_init;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const 
	     typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_Vector>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
get_p_init(LocalOrdinal l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, 
		     std::logic_error,
                     std::endl << 
                     "Error!  twoD_diffusion_problem::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return p_init;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<typename twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,
					     LocalOrdinal,GlobalOrdinal,Node
					     >::Tpetra_CrsMatrix>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
create_W() const
{
  Teuchos::RCP<Tpetra_CrsMatrix> AA = 
    Teuchos::rcp(new Tpetra_CrsMatrix(graph));
  AA->fillComplete();
  return AA;
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void 
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
computeResidual(const Tpetra_Vector& x,
		const Tpetra_Vector& p,
		Tpetra_Vector& f)
{
  // f = A*x - b
  if (log_normal)
    computeA(*lnFunc, p, *A);
  else
    computeA(*klFunc, p, *A);
  A->apply(x,f);
  f.update(-1.0, *b, 1.0);
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void 
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
computeJacobian(const Tpetra_Vector& x,
		const Tpetra_Vector& p,
		Tpetra_CrsMatrix& jac)
{
  if (log_normal)
    computeA(*lnFunc, p, jac);
  else
    computeA(*klFunc, p, jac);
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void 
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
computeResponse(const Tpetra_Vector& x,
		const Tpetra_Vector& p,
		Tpetra_Vector& g)
{
  // g = average of x
  Teuchos::ArrayRCP<Scalar> g_view = g.get1dViewNonConst();
  x.meanValue(g_view());
  g_view[0] *= Scalar(x.getGlobalLength()) / Scalar(mesh.size());
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
template <typename FuncT>
void
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::
computeA(const FuncT& func, const Tpetra_Vector& p, Tpetra_CrsMatrix& jac)
{
  using Teuchos::ArrayView;
  using Teuchos::arrayView;

  jac.resumeFill();
  jac.setAllToScalar(0.0);

  Teuchos::ArrayRCP<const Scalar> p_view = p.get1dView();
  Teuchos::Array<Scalar> rv(p_view());
  size_t NumMyElements = x_map->getLocalNumElements();
  ArrayView<const GlobalOrdinal> MyGlobalElements = 
    x_map->getLocalElementList ();
  MeshScalar h2 = h*h;
  Scalar val;

  for(size_t i=0 ; i<NumMyElements; ++i ) {
      
    // Center
    GlobalOrdinal global_idx = MyGlobalElements[i];
    if (mesh[global_idx].boundary) {
      val = 1.0;
      jac.replaceGlobalValues(global_idx, arrayView(&global_idx,1), 
			      arrayView(&val,1));
    }
    else {
      Scalar a_down = 
	-func(mesh[global_idx].x, mesh[global_idx].y-h/2.0, rv)/h2;
      Scalar a_left = 
	-func(mesh[global_idx].x-h/2.0, mesh[global_idx].y, rv)/h2;
      Scalar a_right = 
	-func(mesh[global_idx].x+h/2.0, mesh[global_idx].y, rv)/h2;
      Scalar a_up = 
	-func(mesh[global_idx].x, mesh[global_idx].y+h/2.0, rv)/h2;
      
      // Center
      val = -(a_down + a_left + a_right + a_up);
      jac.replaceGlobalValues(global_idx, arrayView(&global_idx,1), 
			      arrayView(&val,1));

      // Down
      if (!(eliminate_bcs && mesh[mesh[global_idx].down].boundary))
	jac.replaceGlobalValues(global_idx, 
				arrayView(&mesh[global_idx].down,1), 
				arrayView(&a_down,1));
      
      // Left
      if (!(eliminate_bcs && mesh[mesh[global_idx].left].boundary))
	jac.replaceGlobalValues(global_idx, 
				arrayView(&mesh[global_idx].left,1), 
				arrayView(&a_left,1));
      
      // Right
      if (!(eliminate_bcs && mesh[mesh[global_idx].right].boundary))
	jac.replaceGlobalValues(global_idx, 
				arrayView(&mesh[global_idx].right,1), 
				arrayView(&a_right,1));

      // Up
      if (!(eliminate_bcs && mesh[mesh[global_idx].up].boundary))
	jac.replaceGlobalValues(global_idx, 
				arrayView(&mesh[global_idx].up,1), 
				arrayView(&a_up,1));
    }
  }
  jac.fillComplete();
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::KL_Diffusion_Func::
KL_Diffusion_Func(MeshScalar xyLeft, MeshScalar xyRight, 
		  BasisScalar mean, BasisScalar std_dev, 
		  MeshScalar L, LocalOrdinal num_KL) : point(2)
{
  Teuchos::ParameterList rfParams;
  rfParams.set("Number of KL Terms", num_KL);
  rfParams.set("Mean", mean);
  rfParams.set("Standard Deviation", std_dev);
  LocalOrdinal ndim = 2;
  Teuchos::Array<MeshScalar> domain_upper(ndim), domain_lower(ndim), 
    correlation_length(ndim);
  for (LocalOrdinal i=0; i<ndim; i++) {
    domain_upper[i] = xyRight;
    domain_lower[i] = xyLeft;
    correlation_length[i] = L;
  }
  rfParams.set("Domain Upper Bounds", domain_upper);
  rfParams.set("Domain Lower Bounds", domain_lower);
  rfParams.set("Correlation Lengths", correlation_length);
  
  rf = 
    Teuchos::rcp(new Stokhos::KL::ExponentialRandomField<MeshScalar>(rfParams));
}

template <typename Scalar, typename MeshScalar, typename BasisScalar,
	  typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Scalar 
twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,
		       Node>::KL_Diffusion_Func::
operator() (MeshScalar x, MeshScalar y, const Teuchos::Array<Scalar>& rv) const 
{
  point[0] = x;
  point[1] = y;
  return rf->evaluate(point, rv);
}
