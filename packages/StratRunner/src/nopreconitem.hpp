/*
 * nopreconitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#ifndef NOPRECONITEM_HPP_
#define NOPRECONITEM_HPP_

#include "parameterlisttreeitem.hpp"

/**
 * NoPreconItem is basically a place holder. If the user does not whish to use a preconditioner, this item will be used
 * as the preconditioner.
 *
 * @see MLItem
 * @see IfpackItem
 */
class NoPreconItem: public ParameterListTreeItem{
public:
	enum {Type = UserType +10};
	/**
 	 * Constructs an NoPreconItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the NoPreconItem.
 	 */
	NoPreconItem(ParameterTreeItem *parent);
};
#endif /* NOPRECONNODE_HPP_ */

