/*
 * amesositem.hpp
 *
 *  Created on: Dec 6, 2008
 *      Author: kurtis
 */

#ifndef AMESOSITEM_HPP_
#define AMESOSITEM_HPP_
#include "parameterlisttreeitem.hpp"
/**
 * AmesosItem contains all the input neccessary for using the Amesos Solver in stratimikos.
 *
 * @see AztecOOItem
 * @see BelosItem
 */
class AmesosItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +9};
	/**
 	 * Constructs an AmesosItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the AmesosItem.
 	 */
	AmesosItem(ParameterTreeItem *parent);
private:
	void addParameters();
	void addParameterLists();
};


#endif /* AMESOSITEM_HPP_ */

