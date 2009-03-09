/*
 * mlitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#ifndef MLITEM_HPP_
#define MLITEM_HPP_

#include "parameterlisttreeitem.hpp"

/**
 * MLItem contains all the input neccessary for using the ML Preconditioner in stratimikos.
 *
 * @see IfpackItem
 * @see NoPreconItem
 */
class MLItem : public ParameterListTreeItem{
public:
	enum {Type = UserType + 11};
	/**
 	 * Constructs an MLItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the MLItem.
 	 */
	MLItem(ParameterTreeItem *parent);
private:
	void addParameters();
	void addParameterLists();
};


#endif /* MLNODE_HPP_ */
