/*
 * solveritem.hpp
 *
 *  Created on: Dec 11, 2008
 *      Author: kurtis
 */

#ifndef SOLVERITEM_HPP_
#define SOLVERITEM_HPP_
#include "parameterlisttreeitem.hpp"

/**
 * A special ParameterListTreeItem that allows the user to choose which solver he or she would like to use.
 */
class SolverItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +5};
	/**
	 * Contructs a SolverItem object.
	 *
	 * @param parent The parent ParameterTreeItem.
	 */
	SolverItem(ParameterTreeItem *parent);
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
private:
	ParameterListTreeItem *solverChild;
	QStringList solversList;
};

#endif /* SOLVERITEM_HPP_ */
