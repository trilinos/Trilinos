/*
 * ifpacksettings.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
/*
 * TODO
 * -Figure out the options for sparskit: type
 */
#include "ifpacksettingsitem.hpp"
#include "booltreeitem.hpp"
#include "doubletreeitem.hpp"
#include "inttreeitem.hpp"
#include "stringtreeitem.hpp"

IfpackSettingsItem::IfpackSettingsItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, IfpackSettingsItem::Type)
{
	setText(0,QString("Ifpack Settings"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void IfpackSettingsItem::addParameters(){
	QStringList solverTypeList = QStringList() << "Amesos_Lapack" << "Amesos_Klu" << "Amesos_Umfpack" << "Amesos_Superlu" << "Amesos_Mumps" << "Amesos_Dscpack";
	addChild(new StringTreeItem("amesos: solver type", "Type:", solverTypeList, this));
	addChild(new DoubleTreeItem("fact: absolute threshold", "Threshold:", this));
	addChild(new DoubleTreeItem("fact: drop tolerance", "Tolerance:", this));
	addChild(new DoubleTreeItem("fact: ict level-of-fill", "level-of-fill:", this, 1));
	addChild(new DoubleTreeItem("fact: ilut level-of-fill", "level-of-fill:", this, 1));
	addChild(new IntTreeItem("fact: level-of-fill", "level-of-fill:", this));
	addChild(new DoubleTreeItem("fact: relative threshold", "value:", this, 1));
	addChild(new DoubleTreeItem("fact: relax value", "value:", this));
	addChild(new DoubleTreeItem("fact: sparskit: alph", "value:", this));
	addChild(new DoubleTreeItem("fact: sparskit: droptol", "value:", this));
	addChild(new IntTreeItem("fact: sparskit: lfil", "value:", this));
	addChild(new IntTreeItem("fact: sparskit: mbloc", "value:", this, -1));
	addChild(new DoubleTreeItem("fact: sparskit: permtol", "value:", this, 0.1));
	addChild(new DoubleTreeItem("fact: sparskit: tol", "value:", this));
	//QStringList sparskitTypeList = QStringList() << "ILUT";
	addChild(new IntTreeItem("partitioner: local parts", "value:", this, 1));
	addChild(new IntTreeItem("partitioner: overlap", "value:", this));
	addChild(new IntTreeItem("partitioner: print level", "value:", this));
	QStringList partionerTypeList = QStringList() << "Linear" << "Greedy" << "Metis" << "Equation";
	addChild(new StringTreeItem("partitioner: type", "Type:", partionerTypeList, this));
	addChild(new BoolTreeItem("partitioner: use symmetric graph", this, true));
	addChild(new DoubleTreeItem("relaxation: damping factor:", "value:", this, 1));
	addChild(new DoubleTreeItem("relaxation: min diagonal value", "value:", this, 1));
	addChild(new IntTreeItem("relaxation: sweeps", "value:", this, 1));
	QStringList relaxationTypeList = QStringList() << "Jacobi" << "Gauss-Seidel" << "Symmetric Gauss-Seidel";
	addChild(new StringTreeItem("relaxation: type", "Type:", relaxationTypeList, this));
	addChild(new BoolTreeItem("relaxation: zero starting solution", this));
	QStringList schwarzCombineList = QStringList() << "Zero" << "Add" << "Insert" << "InsertAdd" << "Average" << "AbsMax";
	addChild(new StringTreeItem("schwarz: combine mode", "Mode:", schwarzCombineList, this));
	addChild(new BoolTreeItem("schwarz: compute condest", this, true));
	addChild(new BoolTreeItem("schwarz: filter singletons", this));
	QStringList schwarzReorderingList = QStringList() << "none" << "RCM" << "Metis";
	addChild(new StringTreeItem("schwarz: reordering type", "Type:", schwarzReorderingList, this));
}

