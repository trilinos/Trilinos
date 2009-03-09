/*
 * psuedoblockgmresitem.hpp
 *
 *  Created on: Dec 20, 2008
 *      Author: kurtis
 */

#ifndef PSUEDOBLOCKGMRESITEM_HPP_
#define PSUEDOBLOCKGMRESITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * PsuedoBlockGMRESItem contains all the input neccessary for using the PsuedoBlockGMRES Solver within Belos.
 *
 * @see BelosSolverItem
 */
class PsuedoBlockGMRESItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +22};
	/**
 	 * Constructs an PsuedoBlockGMRESItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the PsuedoBlockGMRESItem.
 	 */
	PsuedoBlockGMRESItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* PSUEDOBLOCKGMRESITEM_HPP_ */
