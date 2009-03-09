/*
 * stratroot.hpp
 *
 *  Created on: Dec 5, 2008
 *      Author: kurtis
 */

#ifndef STRATROOT_HPP_
#define STRATROOT_HPP_
//#include <QTreeWidgetItem>
#include "parameterlisttreeitem.hpp"
class SolverTree;

/**
 * StratRoot is the root item used in a SolverTree. All other ParameterTreeItems should be a children of it.
 */
class StratRoot : public ParameterListTreeItem{
public:
	enum {Type = UserType +1};
	/**
	 * Constructs a StratRoot object.
	 *
	 * @param parent The parent solver tree
	 */
	StratRoot(SolverTree *parent);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
private:
	void addParameters();
	void addParameterLists();
};

#endif /* STRATROOT_HPP_ */
