/*
 * blockcgitem.hpp
 *
 *  Created on: Dec 20, 2008
 *      Author: kurtis
 */

#ifndef BLOCKCGITEM_HPP_
#define BLOCKCGITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * BlockCGItem contains all the input neccessary for using the BlockCG Solver within Belos.
 *
 * @see BelosSolverItem
 */
class BlockCGItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +20};
	/**
 	 * Constructs an BlockCGItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the BlockCGItem.
 	 */
	BlockCGItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* BLOCKCGITEM_HPP_ */
