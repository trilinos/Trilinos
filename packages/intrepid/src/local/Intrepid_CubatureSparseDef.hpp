namespace Intrepid {

/***********************************************************************
**	Helper Functions for Class CubatureSparse
************************************************************************/

//template< class Scalar, int DIM >
//void iterateThroughDimensions(int level, int dims_left, SGNodes< Scalar, DIM > & cubPointsND, Scalar* partial_node, Scalar partial_weight);
template< class Scalar, int DIM >
void iterateThroughDimensions(int level, int dims_left, SGNodes< Scalar, DIM > & cubPointsND, Teuchos::Array<Scalar> & partial_node, Scalar partial_weight); 

int calculateNumPoints(int dim, int level);

//int iterateThroughDimensionsForNumCalc(int dims_left, int level, int levels_left, int level_so_far, int* nodes, int product, bool no_uni_quad);
int iterateThroughDimensionsForNumCalc(int dims_left, int level, int levels_left, int level_so_far, Teuchos::Array<int> & nodes, int product, bool no_uni_quad);

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
		else if(degree <= 19)
			level_ = 6;
		else if(degree <= 23)
			level_ = 7;
		else if(degree <= 27)
			level_ = 8;
		else if(degree <= 31)
			level_ = 9;
		else if(degree <= 35)
			level_ = 10;
		else if(degree <= 39)
			level_ = 11;
		else if(degree <= 43)
			level_ = 12;
		else if(degree <= 47)
			level_ = 13;
		else if(degree <= 51)
			level_ = 14;
		else if(degree <= 55)
			level_ = 15;
		else if(degree <= 59)
			level_ = 16;
	}
	else if(dimension_ == 3)
	{
		if(degree == 1)
			level_ = 1;
		else if(degree <= 3)
			level_ = 2;
		else if(degree <= 5)
			level_ = 3;
		else if(degree <= 9)
			level_ = 4;
		else if(degree <= 13)
			level_ = 5;
		else if(degree <= 17)
			level_ = 6;
		else if(degree <= 21)
			level_ = 7;
		else if(degree <= 25)
			level_ = 8;
		else if(degree <= 29)
			level_ = 9;
		else if(degree <= 33)
			level_ = 10;
		else if(degree <= 37)
			level_ = 11;
		else if(degree <= 41)
			level_ = 12;
		else if(degree <= 45)
			level_ = 13;
		else if(degree <= 49)
			level_ = 14;
		else if(degree <= 53)
			level_ = 15;
		else if(degree <= 57)
			level_ = 16;
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
	//Scalar* dummy_point = {0};
	Teuchos::Array<Scalar> dummy_point(1);
        dummy_point[0] = 0.0;
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

//template< class Scalar, int DIM>
//void iterateThroughDimensions(int level, int dims_left, SGNodes< Scalar, DIM > & cubPointsND, Scalar* partial_node, Scalar partial_weight)
template< class Scalar, int DIM>
void iterateThroughDimensions(int level, int dims_left, SGNodes< Scalar, DIM > & cubPointsND, Teuchos::Array<Scalar> & partial_node, Scalar partial_weight)
{
	int l = level;
	int d = DIM;
	int add_on = d - dims_left;
	int start = dims_left > 1 ? 1 : (int)std::max(1, l);
	int end = l + add_on;

	ECell cellType = CELL_EDGE;

	for(int k_i = start; k_i <= end; k_i++)
	{	/*******************
		**	Slow-Gauss
		********************/
		int order1D = 2*k_i-1;
		/*******************
		**	Fast-Gauss
		********************/
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
			//Scalar* node = new Scalar[d-dims_left+1];
			Teuchos::Array<Scalar> node(d-dims_left+1);
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
				cubPointsND.addNode(&node[0], weight);
			}
		}
	}
}
int calculateNumPoints(int dim, int level)
{
	//int* uninum = new int[level];
	Teuchos::Array<int> uninum(level);
	uninum[0] = 1;
	for(int i = 1; i <= level-1; i++)
	{
		uninum[i] = 2*i;
	}

	int numOfNodes = iterateThroughDimensionsForNumCalc(dim, level, level, 0, uninum, 1, true);
	return numOfNodes;
}

//int iterateThroughDimensionsForNumCalc(int dims_left, int level, int levels_left, int level_so_far, int* nodes, int product, bool no_uni_quad)
int iterateThroughDimensionsForNumCalc(int dims_left, int level, int levels_left, int level_so_far, Teuchos::Array<int> & nodes, int product, bool no_uni_quad)
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

} // end namespace Intrepid
