/*
 * inttreeitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#ifndef INTTREEITEM_HPP_
#define INTTREEITEM_HPP_
#include "parametertreeitem.hpp"

/**
 * IntTreeItem is a ParameterTreeItem that represents a int parameter used in Stratimikos
 */
class IntTreeItem : public ParameterTreeItem{
public:
	enum {Type = UserType+1};
	/**
	 * Constructs a IntTreeItem object.
	 *
	 * @param name The name of the Int parameter (will be displayed to the user).
	 * @param dialogLabel The label that will be displayed next to the input region in the dialog that pops up when
	 * the user attemps to change the value of the parameter.
	 * @param parent The parent ParameterTreeItem.
	 * @param value The default value of the int parameter.
	 * @param validator The validator to be used when validating the input give by the user for this parameter.
	 */
	IntTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent, int value=0, QIntValidator *validator=0);
	/**
	 * Deconstructs a IntTreeItem object.
	 */
	~IntTreeItem();
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void changeValue(int value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
	/*
	 * Gets the minimum value allowed by the validator for this parameter.
	 */
//	int getMin();
	/*
	 * Gets the maximum value allowed by the validator for this parameter.
	 */
//	int getMax();
private:
	QString name, dialogLabel;
	QIntValidator *validator;
};

#endif /* INTTREEITEM_HPP_ */
