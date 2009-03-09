/*
 * aztecooitem.hpp
 *
 *  Created on: Dec 5, 2008
 *      Author: kurtis
 */

#ifndef AZTECOOITEM_HPP_
#define AZTECOOITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * AztecOOItem contains all the input neccessary for using the AztecOO Solver in stratimikos.
 *
 * @see BelosItem
 * @see AmesosItem
 */
class AztecOOItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +8};
	/**
 	 * Constructs an AztecOOItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the AztecOOItem.
 	 */
	AztecOOItem(ParameterTreeItem *parent);
private:
	void addParameters();
	void addParameterLists();
};

#endif /* AZTECOONODE_HPP_ */
