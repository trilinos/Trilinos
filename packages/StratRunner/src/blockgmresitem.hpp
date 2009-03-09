/*
 * blockgmresitem.hpp
 *
 *  Created on: Dec 20, 2008
 *      Author: kurtis
 */

#ifndef BLOCKGMRESITEM_HPP_
#define BLOCKGMRESITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * BlockGMRESItem contains all the input neccessary for using the BlockGMRES Solver within Belos.
 *
 * @see BelosSolverItem
 */
class BlockGMRESItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +21};
	/**
 	 * Constructs an BlockGMRESItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the BlockGMRESItem.
 	 */
	BlockGMRESItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* BLOCKGMRESITEM_HPP_ */

