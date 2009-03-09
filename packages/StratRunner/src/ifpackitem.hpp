/*
 * ifpackitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#ifndef IFPACKITEM_HPP_
#define IFPACKITEM_HPP_

#include "parameterlisttreeitem.hpp"

/**
 * IfpackItem contains all the input neccessary for using the Ifpack Preconditioner in stratimikos.
 *
 * @see MLItem
 * @see NoPreconItem
 */
class IfpackItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +12};
	/**
 	 * Constructs an IfpackItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the IfpackItem.
 	 */
	IfpackItem(ParameterTreeItem *parent);
private:
	void addParameters();
	void addParameterLists();
};


#endif /* IFPACKNODE_HPP_ */
