namespace Intrepid {

/************************************************************************
**	Class Definition for class SGPoint
**	Function: Helper Class with cosntruction of Sparse Grid
*************************************************************************/
template<class Scalar, int D>
class SGPoint{
public:
	Scalar coords[D];

	SGPoint(Scalar p[D]);
	bool const operator==(const SGPoint<Scalar, D> & right);
	bool const operator<(const SGPoint<Scalar, D> & right);
	bool const operator>(const SGPoint<Scalar, D> & right);
	//friend ostream & operator<<(ostream & o, const SGPoint<D> & p);
};

/************************************************************************
**	Class Definition for class SGNode
**	function: Helper Class with constrution of Sparse Grid
************************************************************************/
template<class Scalar, int D>
class SGNodes{
public:
	std::vector< SGPoint<Scalar, D> > nodes;
	std::vector< Scalar > weights;
	bool addNode(Scalar new_node[D], Scalar weight);
	void copyToTeuchos(Teuchos::Array< Point<Scalar> > & cubPoints, Teuchos::Array<Scalar> & cubWeights);
	void copyToTeuchos(Teuchos::Array< Scalar* > & cubPoints, Teuchos::Array<Scalar> & cubWeights);
	int size() {return nodes.size();}
};

template< class Scalar, int DIM >
void iterateThroughDimensions(int level, int dims_left, SGNodes< Scalar, DIM > & cubPointsND, Scalar* partial_node, Scalar partial_weight);

int calculateNumPoints(int dim, int level);
int iterateThroughDimensionsForNumCalc(int dims_left, int level, int levels_left, int level_so_far, int* nodes, int product, bool no_uni_quad);

template< class Scalar, int D >
void getSpecialCub1(SGNodes< Scalar, D > & cubPointsND);


/**************************************************************************
**	Function Definitions for Class CubatureSparse
***************************************************************************/

template <class Scalar, int dimension_>
CubatureSparse<Scalar, dimension_>::CubatureSparse(const int                        degree) :
  degree_(degree) {

	if(dimension_ == 2)
	{
		if(degree == 1)
			level_ = 1;
		else if(degree <= 3)
			level_ = 2;
		else if(degree <= 7)
			level_ = 3;
		else if(degree <= 11)
			level_ = 4;
		else if(degree <= 15)
			level_ = 5;
	}
	else if(dimension_ == 3)
	{
		if(degree <= 1)
			level_ = 1;
		else if(degree <= 3)
			level_ = 2;
		else if(degree <= 5)
			level_ = 3;
		else if(degree <= 9)
			level_ = 4;
		else if(degree <= 13)
			level_ = 5;
	}

	numPoints_ = calculateNumPoints(dimension_,level_);
}



template <class Scalar, int dimension_>
void CubatureSparse<Scalar,dimension_>::getCubature(int &                            numCubPoints,
                                         Teuchos::Array< Point<Scalar> >& cubPoints,
                                         Teuchos::Array<Scalar>&          cubWeights) const {

	CubatureSparse<Scalar,dimension_>::getCubature(cubPoints, cubWeights);
	numCubPoints = numPoints_;
} // end getCubature



template <class Scalar, int dimension_>
void CubatureSparse<Scalar,dimension_>::getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                                         Teuchos::Array<Scalar>&          cubWeights) const{
	Scalar* dummy_point = {0};
	Scalar dummy_weight = 1.0;
	SGNodes<Scalar, dimension_> grid;

	iterateThroughDimensions(level_, dimension_, grid, dummy_point, dummy_weight);
	grid.copyToTeuchos(cubPoints, cubWeights);
} // end getCubature

template<class Scalar, int dimension_>
int CubatureSparse<Scalar, dimension_>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints

template <class Scalar, int dimension_>
ECell CubatureSparse<Scalar,dimension_>::getCellType() const {

	ECell celltype = CELL_HEX;  

	if(dimension_ == 2)
		celltype = CELL_QUAD;

	return celltype;
} // end getCellType



template <class Scalar,int dimension_>
int CubatureSparse<Scalar,dimension_>::getAccuracy() const {
  return 0;
} //end getAccuracy

template <class Scalar,int dimension_>
void CubatureSparse<Scalar,dimension_>::getSpecialCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const{
	SGNodes<Scalar, 2> grid;

	getSpecialCub1(grid);
	grid.copyToTeuchos(cubPoints, cubWeights);
}


/**************************************************************************
**	Function Definitions for Class SGPoint
***************************************************************************/
template<class Scalar, int D>
SGPoint<Scalar, D>::SGPoint(Scalar p[D])
{
	for(int i = 0; i < D; i++)
	{
		coords[i] = p[i];
	}
}

template<class Scalar, int D>
bool const SGPoint<Scalar, D>::operator==(const SGPoint<Scalar, D> & right)
{
	bool equal = true;

	for(int i = 0; i < D; i++)
	{
		if(coords[i] != right.coords[i])
			return false;
	}

	return equal;
}

template<class Scalar, int D>
bool const SGPoint<Scalar, D>::operator<(const SGPoint<Scalar, D> & right)
{
	for(int i = 0; i < D; i++)
	{
		if(coords[i] < right.coords[i])
			return true;
		else if(coords[i] > right.coords[i])
			return false;
	}
	
	return false;
}

template<class Scalar, int D>
bool const SGPoint<Scalar, D>::operator>(const SGPoint<Scalar, D> & right)
{
	if(this < right || this == right)
		return false;

	return true;
}

template<class Scalar, int D>
std::ostream & operator<<(std::ostream & o, SGPoint<Scalar, D> & p)
{
	o << "(";
	for(int i = 0; i<D;i++)
		o<< p.coords[i] << " ";
	o << ")";
	return o;
}


/**************************************************************************
**	Function Definitions for Class SGNode
***************************************************************************/

template<class Scalar, int D>
bool SGNodes< Scalar, D >::addNode(Scalar new_node[D], Scalar weight)
{
	SGPoint<Scalar, D> new_point(new_node);
	bool new_and_added = true;

	if(nodes.size() == 0)
	{
		nodes.push_back(new_point);
		weights.push_back(weight);
	}
	else
	{		
		int left = -1;
		int right = (int)nodes.size();
		int mid_node = (int)ceil(nodes.size()/2.0)-1;

		bool iterate_continue = true;

		while(iterate_continue)
		{
			if(new_point == nodes[mid_node]){
				weights[mid_node] += weight;
				iterate_continue = false;
				new_and_added = false;	
			}
			else if(new_point < nodes[mid_node]){
				if(right - left <= 2)
				{
					//insert the node into the vector to the left of mid_node
					nodes.insert(nodes.begin()+mid_node, new_point);
					weights.insert(weights.begin()+mid_node,weight);
					iterate_continue = false;
				}
				else 
				{
					right = mid_node;
					mid_node += (int)ceil((left-mid_node)/2.0);
				}
			}
			else{ //new_point > nodes[mid_node];

				if(mid_node == (int)nodes.size()-1)
				{
					nodes.push_back(new_point);
					weights.push_back(weight);
					iterate_continue = false;
				}
				else if(right - left <= 2)
				{
					//insert the node into the vector to the right of mid_node
					nodes.insert(nodes.begin()+mid_node+1, new_point);
					weights.insert(weights.begin()+mid_node+1,weight);
					iterate_continue = false;
				}
				else 
				{
					left = mid_node;
					mid_node += (int)ceil((right-mid_node)/2.0);
				}
			}
		}
	}

	return new_and_added;
}

template< class Scalar, int D >
void SGNodes<Scalar, D>::copyToTeuchos(Teuchos::Array< Point<Scalar> > & cubPoints, Teuchos::Array<Scalar> & cubWeights)
{
	int numPoints = size();

	Point<Scalar> tempPoint(D);
	cubPoints.assign(numPoints,tempPoint);
	cubWeights.assign(numPoints, 0.0);

	for(int i = 0; i < numPoints; i++)
	{
		cubPoints[i].setCoordinates(nodes[i].coords, D);
		cubWeights[i] = weights[i];
	}
}

template< class Scalar, int D >
void SGNodes<Scalar,D>::copyToTeuchos(Teuchos::Array< Scalar* > & cubPoints, Teuchos::Array<Scalar> & cubWeights)
{
	int numPoints = size();

	Scalar tempPoint[D];
	cubPoints.assign(numPoints,tempPoint);
	cubWeights.assign(numPoints, 0.0);

	for(int i = 0; i < numPoints; i++)
	{
		cubPoints[i] = nodes[i].coords;
		cubWeights[i] = weights[i];
	}
}


/************************************************************************
**Function Definition for iterateThroughDimensions() 
**		and its helper functions
*************************************************************************/
int factorial(int num)
{
	int answer = 1;
	if(num >= 1)
	{
		while(num > 0)
		{
			answer = answer*num;
			num--;
		}
	}
	else if(num == 0)
		answer = 1;
	else
		answer = -1;

	return answer; 
}

double combination(int top, int bot)
{
	double answer = factorial(top)/(factorial(bot) * factorial(top-bot));
	return answer;
}

template< class Scalar, int DIM>
void iterateThroughDimensions(int level, int dims_left, SGNodes< Scalar, DIM > & cubPointsND, Scalar* partial_node, Scalar partial_weight)
{
	int l = level;
	int d = DIM;
	int add_on = d - dims_left;
	int start = dims_left > 1 ? 1 : (int)std::max(1, l);
	int end = l + add_on;

	ECell cellType = CELL_EDGE;

	for(int k_i = start; k_i <= end; k_i++)
	{
		int order1D = 2*k_i-1;	//line that could be more generalized - depending on speed
		//int order1D = (int)pow(2,k_i) - 1;
		Teuchos::Array< Point<Scalar> > cubPoints1D;
		Teuchos::Array<Scalar> cubWeights1D;
		Point<Scalar> tempPoint1D(1);
		cubPoints1D.assign(order1D,tempPoint1D);
		cubWeights1D.assign(order1D, 0.0);
		int cubDegree1D = 2*order1D - 1;
		CubatureDirect<Scalar> Cub1D(cellType, cubDegree1D);
		Cub1D.getCubature(order1D, cubPoints1D, cubWeights1D);

		for(int node1D = 0; node1D < order1D; node1D++)
		{
			Scalar* node = new Scalar[d-dims_left+1];
			Teuchos::Array<Scalar> temp_array = cubPoints1D[node1D].getCoordinates();
			node[d-dims_left] = temp_array[0];
			for(int m = 0; m < d-dims_left; m++)
				node[m] = partial_node[m];

			Scalar weight = cubWeights1D[node1D]*partial_weight;

			if(dims_left != 1)
			{
				iterateThroughDimensions(l - k_i, dims_left-1, cubPointsND, node, weight);
			}
			else
			{
				weight = pow(-1, end - k_i)*combination(d-1, k_i - l)*weight;
				cubPointsND.addNode(node, weight);
			}
		}
	}
}
int calculateNumPoints(int dim, int level)
{
	int* uninum = new int[level];
	uninum[0] = 1;
	for(int i = 1; i <= level-1; i++)
	{
		uninum[i] = 2*i;
	}

	int numOfNodes = iterateThroughDimensionsForNumCalc(dim, level, level, 0, uninum, 1, true);
	return numOfNodes;
}

int iterateThroughDimensionsForNumCalc(int dims_left, int level, int levels_left, int level_so_far, int* nodes, int product, bool no_uni_quad)
{
	int numNodes = 0;
	for(int j = 1; j <= levels_left; j++)
	{
		bool temp_bool = no_uni_quad;
		int temp_knots = nodes[j-1]*product;
		int temp_lsf = level_so_far + j;

		if(j==1)
			temp_bool = false;

		if(dims_left == 1)
		{
			if(temp_lsf < level && temp_bool == true)
				numNodes += 0;
			else
			{
				numNodes += temp_knots;
			}
		}
		else
		{
			numNodes += iterateThroughDimensionsForNumCalc(dims_left-1,level, levels_left-j+1, temp_lsf, nodes, temp_knots, temp_bool);
		}
	}
	return numNodes;
}

template< class Scalar, int D>
void getSpecialCub1(SGNodes< Scalar, D > & cubPointsND)
{
	ECell cellType = CELL_EDGE;
	Point<Scalar> tempPoint1D(1);
	int level1, level2;
	int order1, order2;
	int degree1, degree2;
	
	//Step 1
	level1 = 3;
	level2 = 2;
	order1 = 2*level1-1;
	order2 = 2*level2-1;
	degree1 = 2*order1-1;
	degree2 = 2*order2-1;

	Teuchos::Array< Point<Scalar> > cubPoints1A;
	Teuchos::Array<Scalar> cubWeights1A;
	cubPoints1A.assign(order1,tempPoint1D);
	cubWeights1A.assign(order1, 0.0);
	CubatureDirect<Scalar> Cub1A(cellType, degree1);
	Cub1A.getCubature(order1, cubPoints1A, cubWeights1A);

	Teuchos::Array< Point<Scalar> > cubPoints1B;
	Teuchos::Array<Scalar> cubWeights1B;
	cubPoints1B.assign(order2,tempPoint1D);
	cubWeights1B.assign(order2, 0.0);
	CubatureDirect<Scalar> Cub1B(cellType, degree2);
	Cub1B.getCubature(order2, cubPoints1B, cubWeights1B);

	for(int i = 0; i < order1; i++)
	{
		for(int j = 0; j < order2; j++)
		{
			Scalar node[2];
			Teuchos::Array<Scalar> temp_array = cubPoints1A[i].getCoordinates();
			node[0] = temp_array[0];
			temp_array = cubPoints1B[j].getCoordinates();
			node[1] = temp_array[0];

			Scalar weight = cubWeights1A[i]*cubWeights1B[j];

			cubPointsND.addNode(node, weight);
		}
	}

	//Step 2
	level1 = 2;
	level2 = 3;
	order1 = 2*level1-1;
	order2 = 2*level2-1;
	degree1 = 2*order1-1;
	degree2 = 2*order2-1;

	Teuchos::Array< Point<Scalar> > cubPoints1C;
	Teuchos::Array<Scalar> cubWeights1C;
	cubPoints1C.assign(order1,tempPoint1D);
	cubWeights1C.assign(order1, 0.0);
	CubatureDirect<Scalar> Cub1C(cellType, degree1);
	Cub1C.getCubature(order1, cubPoints1C, cubWeights1C);

	Teuchos::Array< Point<Scalar> > cubPoints1D;
	Teuchos::Array<Scalar> cubWeights1D;
	cubPoints1D.assign(order2,tempPoint1D);
	cubWeights1D.assign(order2, 0.0);
	CubatureDirect<Scalar> Cub1D(cellType, degree2);
	Cub1D.getCubature(order2, cubPoints1D, cubWeights1D);

	for(int i = 0; i < order1; i++)
	{
		for(int j = 0; j < order2; j++)
		{
			Scalar node[2];
			Teuchos::Array<Scalar> temp_array = cubPoints1C[i].getCoordinates();
			node[0] = temp_array[0];
			temp_array = cubPoints1D[j].getCoordinates();
			node[1] = temp_array[0];

			Scalar weight = cubWeights1C[i]*cubWeights1D[j];

			cubPointsND.addNode(node, weight);
		}
	}

	//Step 3
	level1 = 2;
	level2 = 2;
	order1 = 2*level1-1;
	order2 = 2*level2-1;
	degree1 = 2*order1-1;
	degree2 = 2*order2-1;

	Teuchos::Array< Point<Scalar> > cubPoints1E;
	Teuchos::Array<Scalar> cubWeights1E;
	cubPoints1E.assign(order1,tempPoint1D);
	cubWeights1E.assign(order1, 0.0);
	CubatureDirect<Scalar> Cub1E(cellType, degree1);
	Cub1E.getCubature(order1, cubPoints1E, cubWeights1E);

	Teuchos::Array< Point<Scalar> > cubPoints1F;
	Teuchos::Array<Scalar> cubWeights1F;
	cubPoints1F.assign(order2,tempPoint1D);
	cubWeights1F.assign(order2, 0.0);
	CubatureDirect<Scalar> Cub1F(cellType, degree2);
	Cub1F.getCubature(order2, cubPoints1F, cubWeights1F);

	for(int i = 0; i < order1; i++)
	{
		for(int j = 0; j < order2; j++)
		{
			Scalar node[2];
			Teuchos::Array<Scalar> temp_array = cubPoints1E[i].getCoordinates();
			node[0] = temp_array[0];
			temp_array = cubPoints1F[j].getCoordinates();
			node[1] = temp_array[0];

			Scalar weight = -cubWeights1E[i]*cubWeights1F[j];

			cubPointsND.addNode(node, weight);
		}
	}
}

/*template< class Scalar, int D>
void addProduct(int level_part[D], int dims_left, SGNodes< Scalar, D > & cubPointsND, Scalar* partial_node, Scalar partial_weight)
{
	int order[D];
	int cub_degree[D];
	ECell cellType = CELL_EDGE;
	Teuchos::Array< Point<Scalar> > cubPoints1D[D];
	Teuchos::Array<Scalar> cubWeights1D[D];
	Point<Scalar> tempPoint1D(1);

	int totalnum = 1;

	for(int i = 0; i < D; i++)
	{
		totalnum *= level_part[i];

		order[i] = 2*level_part[i] - 1;
		cub_degree[i] = 2*order[i] - 1;

		cubPoints1D[i].assign(order,tempPoint1D);
		cubWeights1D[i].assign(order, 0.0);
		CubatureDirect<Scalar> Cub1D(cellType, cub_degree[i]);
		Cub1D.getCubature(order[i], cubPoints1D[i], cubWeights1D[i]);
	}

	
	Scalar [][] point = new Scalar[D][totalnum];
	Scalar weight[D];

	for(int i = 0; i < totalnum; i++)
	{
		for(int d = 0; d < D; d++)
		{

		}
	}
}*/

} // end namespace Intrepid
