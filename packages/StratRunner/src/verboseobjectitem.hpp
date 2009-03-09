/*
 * verboseobjectitem.hpp
 *
 *  Created on: Dec 19, 2008
 *      Author: kurtis
 */

#ifndef VERBOSEOBJECTITEM_HPP_
#define VERBOSEOBJECTITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * VerboseObjectItem contains all the input neccessary for using a VerboseObject parameter in stratimikos.
 */
class VerboseObjectItem : public ParameterListTreeItem{
public:
	enum {Type = UserType+17};
	/**
 	 * Constructs an VerboseObjectItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the VerboseObjectItem.
 	 */
	VerboseObjectItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* VERBOSEOBJECTITEM_HPP_ */
