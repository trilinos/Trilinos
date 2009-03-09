/*
 * amesositem.cpp
 *
 *  Created on: Dec 6, 2008
 *      Author: kurtis
 */
#include "amesositem.hpp"
#include "stringtreeitem.hpp"
#include "booltreeitem.hpp"
#include "amesossettingsitem.hpp"
#include "verboseobjectitem.hpp"
#include "amesossettingsitem.hpp"

AmesosItem::AmesosItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, AmesosItem::Type)
{
	setText(0,QString("Amesos"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void AmesosItem::addParameters(){
	QStringList refactChoices = QStringList() << "RepivotOnRefactorization" << "NoPivotOnRefactorization";
	addChild(new StringTreeItem("Refactorization Policy", "Policy:", refactChoices, this));

	QStringList typeChoices = QStringList() << "Lapack" << "Klu" << "Umfpack" << "Superlu" << "Superludist" << "Taucs" << "Pardiso" << "Pastix" << "Paraklete" << "Mumps"  << "Scalapack" << "Dscpack";
	addChild(new StringTreeItem("Solver Type", "Type:", typeChoices, this));

	addChild(new BoolTreeItem("Throw On Preconditioner Input", this));
}

void AmesosItem::addParameterLists(){
	addChild(new AmesosSettingsItem(this));
	addChild(new VerboseObjectItem(this));
}

