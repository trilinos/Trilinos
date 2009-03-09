/*
 * belositem.hpp
 *
 *  Created on: Dec 5, 2008
 *      Author: kurtis
 */

#ifndef BELOSITEM_HPP_
#define BELOSITEM_HPP_

#include "parameterlisttreeitem.hpp"
/**
 * BelosItem contains all the input neccessary for using the Amesos Solver in stratimikos.
 *
 * @see AztecOOItem
 * @see AmesosItem
 */
class BelosItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +7};
	/**
 	 * Constructs an BelosItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the BelosItem.
 	 */
	BelosItem(ParameterTreeItem *parent);
private:
	void addParameterLists();
};

#endif /* BELOSNODE_HPP_ */
