/*
 * booltreeitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

/*
 * TODO
 * -instead of using value, use the setData fucntion to store the value like it should be, within the QTreeWidgetItem
 */
#ifndef BOOLTREEITEM_HPP_
#define BOOLTREEITEM_HPP_
#include "parametertreeitem.hpp"

/**
 * BoolTreeItem is a ParameterTreeItem that represents a boolean parameter used in Stratimikos
 */
class BoolTreeItem : public ParameterTreeItem{
public:
	enum {Type = UserType+1};
	/**
	 * Constructs a BoolTreeItem object.
	 *
	 * @param name The name of the boolean parameter (will be displayed to the user).
	 * @param parent The parent ParameterTreeItem.
	 * @param value The default value of the boolean parameter.
	 */
	BoolTreeItem(QString name, ParameterTreeItem *parent, bool value=false);
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void changeValue(bool value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
private:
	QString name, dialogLabel;
	//bool value;
};


#endif /* BOOLTREEITEM_HPP_ */
