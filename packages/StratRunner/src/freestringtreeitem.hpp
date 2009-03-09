/*  freestringtreeitem.hpp
 *
 *  Created on: Feb 10, 2009
 *      Author: kurtis
 */

#ifndef FREESTRINGTREEITEM_HPP_
#define FREESTRINGTREEITEM_HPP_
#include "parametertreeitem.hpp"

/**
 * FreeStringTreeItem is a ParameterTreeItem that represents a string item that can take any sort of input used in Stratimikos
 */
class FreeStringTreeItem : public ParameterTreeItem{
public:
	enum {Type = UserType+23};
	/**
	 * Constructs a FreeStringTreeItem object.
	 *
	 * @param name The name of the FreeString parameter (will be displayed to the user).
	 * @param dialogLabel The label that will be displayed next to the input region in the dialog that pops up when
	 * the user attemps to change the value of the parameter.
	 * @param parent The parent ParameterTreeItem.
	 * @param default The default value of the int parameter.
	 */
	FreeStringTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent, QString defaultValue=QString(""));
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
};

#endif /* FREESTRINGTREEITEM_HPP_ */
