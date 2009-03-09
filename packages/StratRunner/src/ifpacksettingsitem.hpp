/*
 * ifpacksettingsitem.hpp
 *
 *  Created on: Dec 22, 2008
 *      Author: kurtis
 */

#ifndef IFPACKSETTINGSITEM_HPP_
#define IFPACKSETTINGSITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * IfpackSettingsItem contains all the inputs neccesary for specifing the settings to be used with an Ifpack preconditioner.
 */
class IfpackSettingsItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +15};
	/**
	 * Constructs an IfpackSettingsItem object.
	 *
	 * @param parent The parent ParameterTreeItem of the IfpackSettingsItem
	 */
	IfpackSettingsItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* IFPACKSETTINGSITEM_HPP_ */
