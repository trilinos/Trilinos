/*
 * preconitem.hpp
 *
 *  Created on: Dec 11, 2008
 *      Author: kurtis
 */

#ifndef PRECONITEM_HPP_
#define PRECONITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * A special ParameterListTreeItem that allows the user to choose which preconditioner he or she would like to use.
 */
class PreconItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +6};
	/**
	 * Contructs a PreconItem object.
	 *
	 * @param parent The parent ParameterTreeItem.
	 */
	PreconItem(ParameterTreeItem *parent);
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
private:
	ParameterListTreeItem *preconChild;
	QStringList preconsList;
};

#endif /* PRECONITEM_HPP_ */
