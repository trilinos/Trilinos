/*
 * aztecoosettingsitem.hpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */

#ifndef AZTECOOSETTINGSITEM_HPP_
#define AZTECOOSETTINGSITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * AztecOOSettingsItem contains all the inputs neccesary for specifing the settings to be used with an AztecOO Solver
 */
class AztecOOSettingsItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +13};
	/**
	 * Constructs an AztecOOSettingsItem object.
	 *
	 * @param parent The parent ParameterTreeItem of the AztecOOSettingsItem
	 */
	AztecOOSettingsItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* AZTECOOSETTINGSITEM_HPP_ */
