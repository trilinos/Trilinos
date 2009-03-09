/*
 * belossolveritem.hpp
 *
 *  Created on: Dec 20, 2008
 *      Author: kurtis
 */

#ifndef BELOSSOLVERITEM_HPP_
#define BELOSSOLVERITEM_HPP_
#include "parameterlisttreeitem.hpp"


/**
 * Allows the user to choose with solver to use for the Belos Solver.
 *
 * @see BlockGMRESItem
 * @see BlockCGItem
 * @see PsuedoBlockGMRESItem
 */
class BelosSolverItem : public ParameterListTreeItem{
public:
	enum {Type = UserType +19};
	/**
 	 * Constructs an BelosSolverItem object.
 	 *
 	 * @param parent The parent ParameterTreeItem of the BelosSolverItem.
 	 */
	BelosSolverItem(ParameterTreeItem *parent);
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

#endif /* BELOSSOLVERITEM_HPP_ */
