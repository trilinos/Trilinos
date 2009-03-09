/*
 * amesossettingsitem.hpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */

#ifndef AMESOSSETTINGSITEM_HPP_
#define AMESOSSETTINGSITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * AmesosSettingsItem contains all the inputs neccesary for specifing the settings to be used with an Amesos Solver
 */
class AmesosSettingsItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +14};
	/**
	 * Constructs an AmesosSettingsItem object.
	 *
	 * @param parent The parent ParameterTreeItem of the AmesosSettingsItem
	 */
	AmesosSettingsItem(ParameterTreeItem *parent);
private:
	void addParameters();
	void addParameterLists();
};

#endif /* AMESOSSETTINGSITEM_HPP_ */
