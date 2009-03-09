/*
 * stringtreeitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#ifndef STRINGTREEITEM_HPP_
#define STRINGTREEITEM_HPP_
#include "parametertreeitem.hpp"

/*
 * TODO
 * -instead of using value, use the setData fucntion to store the value like it should be, within the QTreeWidgetItem
 */

/**
 * StringTreeItem is a ParameterTreeItem that represents a string parameter used in Stratimikos
 */
class StringTreeItem : public ParameterTreeItem{
public:
	enum {Type = UserType+4};
	/**
	 * Constructs a BoolTreeItem object.
	 *
	 * @param name The name of the boolean parameter (will be displayed to the user).
	 * @param dialogLabel The label that will be displayed next to the input region in the dialog that pops up when the user attemps to change the value of the parameter.
	 * @param validChoices A list of valid choices for the parameter.
	 * @param parent The parent ParameterTreeItem.
	 */
	StringTreeItem(QString name, QString dialogLabel, QStringList validChoices, ParameterTreeItem *parent);
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void changeValue(QString value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
private:
	QString name, dialogLabel;
	//QString value;
	QStringList validChoices;
};

#endif /* STRINGTREEITEM_HPP_ */
