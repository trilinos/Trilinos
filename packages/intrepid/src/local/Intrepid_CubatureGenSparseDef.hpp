#include "Intrepid_CubatureSparseHelper.hpp"

namespace Intrepid {
/**************************************************************************
**	Definitions for Helper Functions
***************************************************************************/

template<class Scalar>
Scalar Sum(Scalar* list, int first, int last);

/**************************************************************************
**	Function Definitions for Class CubatureSparse
***************************************************************************/

template <class Scalar, int dimension_>
CubatureGenSparse<Scalar, dimension_>::CubatureGenSparse(const int                        degree) :
  degree_(degree) {
SGNodes<int, dimension_> list;
	SGNodes<int,dimension_> bigger_rules;

	bool continue_making_first_list = true;
	bool more_bigger_rules = true;

	ECell cellType = CELL_EDGE;

	int poly_exp[dimension_];
	int level[dimension_];
	int temp_big_rule[dimension_];
	
	for(int i = 0; i<dimension_; i++){
		poly_exp[i] = 0;
		temp_big_rule[i] = 0;
	}

	while(continue_making_first_list){
		for(int i = 0; i < dimension_; i++)
		{
			int max_exp = 0;
			if(i == 0)
				max_exp = std::max(degree_,1) - Sum(poly_exp,1,dimension_-1);
			else if(i == dimension_ -1)
				max_exp = std::max(degree_,1) - Sum(poly_exp,0,dimension_-2);
			else
				max_exp = std::max(degree_,1) - Sum(poly_exp,0,dimension_-1) + poly_exp[i];

			if(poly_exp[i] < max_exp)
			{
				poly_exp[i]++;
				break;
			}
			else
			{
				if(i == dimension_-1)
					continue_making_first_list = false;
				else
					poly_exp[i] = 0;
					
			}
		}

		if(continue_making_first_list)
		{
			for(int j = 0; j < dimension_;j++)
			{
				/*******************
				**	Slow-Gauss
				********************/
				level[j] = (int)std::ceil((((double)poly_exp[j])+3.0)/4.0);
				/*******************
				**	Fast-Gauss
				********************/
				//level[j] = intstd::ceil(std::log(poly_exp[j]+3)/std::log(2) - 1);
			}
			list.addNode(level,1);
			
		}
	}

	/*std::cout << "List:\n";

	for(int i = 0; i<list.size(); i++)
	{
		std::cout << list.nodes[i] << "\n";
	}*/

	while(more_bigger_rules)
	{
		bigger_rules.addNode(temp_big_rule,1);

		for(int i = 0; i < dimension_; i++)
		{
			if(temp_big_rule[i] == 0){
				temp_big_rule[i] = 1;
				break;
			}
			else{
				if(i == dimension_-1)
					more_bigger_rules = false;
				else
					temp_big_rule[i] = 0;
			}
		}
	}	

	for(int x = 0; x < list.size(); x++){
		for(int y = 0; y < bigger_rules.size(); y++)
		{	
			SGPoint<int, dimension_> next_rule;
			for(int t = 0; t < dimension_; t++)
				next_rule.coords[t] = list.nodes[x].coords[t] + bigger_rules.nodes[y].coords[t];

			bool is_in_set = false;
			for(int z = 0; z < list.size(); z++)
			{
				if(next_rule == list.nodes[z]){
					is_in_set = true;
					break;
				}
			}

			if(is_in_set)
			{
				int big_sum[dimension_];
				for(int i = 0; i < dimension_; i++)
					big_sum[i] = bigger_rules.nodes[y].coords[i];
				double coeff = std::pow(-1.0, Sum(big_sum, 0, dimension_-1));
				
				double point[dimension_];
				int point_record[dimension_];

				for(int j = 0; j<dimension_; j++)
					point_record[j] = 1;

				bool more_points = true;

				while(more_points)
				{
					double weight = 1.0;
				
					for(int w = 0; w < dimension_; w++){
						/*******************
						**	Slow-Gauss
						********************/
						int order1D = 2*list.nodes[x].coords[w]-1;
						/*******************
						**	Fast-Gauss
						********************/
						//int order1D = (int)std::pow(2.0,next_rule.coords[w]) - 1;
						Teuchos::Array< Point<Scalar> > cubPoints1D;
						Teuchos::Array<Scalar> cubWeights1D;
						Point<Scalar> tempPoint1D(1);
						cubPoints1D.assign(order1D,tempPoint1D);
						cubWeights1D.assign(order1D, 0.0);
						int cubDegree1D = 2*order1D - 1;
						CubatureDirect<Scalar> Cub1D(cellType, cubDegree1D);
						Cub1D.getCubature(order1D, cubPoints1D, cubWeights1D);

						Teuchos::Array<Scalar> temp_array = cubPoints1D[point_record[w]-1].getCoordinates();
						point[w] = temp_array[0];
						weight = weight * cubWeights1D[point_record[w]-1];
					}			
					weight = weight*coeff;
					grid.addNode(point, weight);

					for(int v = 0; v < dimension_; v++)
					{
						if(point_record[v] < 2*list.nodes[x].coords[v]-1){
							(point_record[v])++;
							break;
						}
						else{
							if(v == dimension_-1)
								more_points = false;
							else
								point_record[v] = 1;
						}
					}
				}
			}
		}
	}

	numPoints_ = grid.size();
}

template <class Scalar, int dimension_>
void CubatureGenSparse<Scalar,dimension_>::getCubature(int &                            numCubPoints,
				Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const {

	CubatureGenSparse<Scalar,dimension_>::getCubature(cubPoints, cubWeights);
	numCubPoints = numPoints_;
} // end getCubature

template <class Scalar, int dimension_>
void CubatureGenSparse<Scalar,dimension_>::getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                           Teuchos::Array<Scalar>&          cubWeights) const{
	grid.copyToTeuchos(cubPoints, cubWeights);
} // end getCubature

/**************************************************************************
**	Function Implementations for Helper Functions
***************************************************************************/

template<class Scalar>
Scalar Sum(Scalar* list, int first, int last)
{
	Scalar sum = 0;
	for(int i = first; i <= last; i++)
		sum += list[i];
	return sum;
}

} // end namespace Intrepid
