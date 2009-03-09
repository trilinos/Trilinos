/*
 * doubletreeitem.hpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#ifndef DOUBLETREEITEM_HPP_
#define DOUBLETREEITEM_HPP_
#include "parametertreeitem.hpp"

/**
 * DoubleTreeItem is a ParameterTreeItem that represents a double parameter used in Stratimikos
 */
class DoubleTreeItem : public ParameterTreeItem{
public:
	enum {Type = UserType+2};
	/**
	 * Constructs a DoubleTreeItem object.
	 *
	 * @param name The name of the Double parameter (will be displayed to the user).
	 * @param dialogLabel The label that will be displayed next to the input region in the dialog that pops up when the user attemps to change the value of the parameter.
	 * @param parent The parent ParameterTreeItem.
	 * @param value The default value of the double parameter.
	 * @param validator The validator to be used when validating the input give by the user for this parameter.
	 */
	DoubleTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent, double value=0, QDoubleValidator *validator=0);
	/**
	 * Deconstructs a DoubleTreeItem object.
	 */
	~DoubleTreeItem();
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void changeValue(double value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
	/**
	 * Gets the current value of this parameter.
	 */
	//double getValue();
	/**
	 * Gets the minimum value allowed by the validator for this parameter.
	 */
	//double getMin();
	/**
	 * Gets the maximum value allowed by the validator for this parameter.
	 */
	//double getMax();
private:
	QString name, dialogLabel;
	//double value;
	QDoubleValidator *validator;
};

#endif /* DOUBLETREEITEM_HPP_ */
