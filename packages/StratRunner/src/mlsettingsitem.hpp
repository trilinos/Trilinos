/*
 * mlsettingsitem.hpp
 *
 *  Created on: Dec 22, 2008
 *      Author: kurtis
 */

#ifndef MLSETTINGSITEM_HPP_
#define MLSETTINGSITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * MLSettingsItem contains all the inputs neccesary for specifing the settings to be used with an ML preconditioner.
 */
class MLSettingsItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +15};
	/**
	 * Constructs an MLSettingsItem object.
	 *
	 * @param parent The parent ParameterTreeItem of the MLSettingsItem
	 */
	MLSettingsItem(ParameterTreeItem *parent);
private:
	void addParameters();
};

#endif /* MLSETTINGSITEM_HPP_ */
