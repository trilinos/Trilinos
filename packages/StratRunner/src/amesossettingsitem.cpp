/*
 * amesossettings.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
#include "amesossettingsitem.hpp"
#include "booltreeitem.hpp"
#include "doubletreeitem.hpp"
#include "inttreeitem.hpp"
#include "stringtreeitem.hpp"

AmesosSettingsItem::AmesosSettingsItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, AmesosSettingsItem::Type)
{
	setText(0,QString("Amesos Settings"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void AmesosSettingsItem::addParameters(){
	addChild(new BoolTreeItem("AddZeroToDiag", this));
	addChild(new BoolTreeItem("ComputeTrueResidual", this));
	addChild(new BoolTreeItem("ComputeVectorNorms", this));
	addChild(new BoolTreeItem("NoDestroy", this));
	addChild(new BoolTreeItem("PrintTiming", this));
	addChild(new BoolTreeItem("Redistribute", this));
	addChild(new BoolTreeItem("Refactorize", this));
	addChild(new BoolTreeItem("TrustMe", this));

	addChild(new DoubleTreeItem("AddToDiag", "Add:", this));
	addChild(new DoubleTreeItem("RcondThreshold", "Threshold:", this));

	addChild(new IntTreeItem("DebugLevel", "Level:", this));
	addChild(new IntTreeItem("MaxProcs", "Maximum:", this));
	addChild(new IntTreeItem("OutputLevel", "Level:", this));
	addChild(new IntTreeItem("Reindex", "Reindex:", this));
	addChild(new IntTreeItem("ScaleMethod", "Method:", this));

	QStringList properties = QStringList() << "General Unsymetric Matrix" << "General Symetric Matrix" << "SPD";
	addChild(new StringTreeItem("MatrixProperty", "Property:", properties, this));


}

void AmesosSettingsItem::addParameterLists(){
	QList<ParameterTreeItem*> lapackItems;
	lapackItems.append(new BoolTreeItem("Equilibrate", this));
	ParameterListTreeItem *lapackList = new ParameterListTreeItem("Lapack", this, lapackItems);
	addChild(lapackList);

	QList<ParameterTreeItem*> mumpsItems;
	mumpsItems.append(new BoolTreeItem("Equilibrate", this, true));
	mumpsItems.append(new DoubleTreeItem("ColScaling", "Scailing:", this));
	mumpsItems.append(new DoubleTreeItem("RowScaling", "Scailing:", this));
	ParameterListTreeItem *mumpsList = new ParameterListTreeItem("Mumps", this, mumpsItems);
	addChild(mumpsList);

	QList<ParameterTreeItem*> pardisoItems;
	pardisoItems.append(new IntTreeItem("IFPARM(1)", "IFPARM(1):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(2)", "IFPARM(2):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(3)", "IFPARM(3):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(4)", "IFPARM(4):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(8)", "IFPARM(8):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(10)", "IFPARM(10):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(11)", "IFPARM(11):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(18)", "IFPARM(18):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(19)", "IFPARM(19):", this));
	pardisoItems.append(new IntTreeItem("IFPARM(21)", "IFPARM(21):", this));
	pardisoItems.append(new IntTreeItem("MSGLVL", "MSGLVL:", this));
	ParameterListTreeItem *pardisoList = new ParameterListTreeItem("Pardiso", this, pardisoItems);
	addChild(pardisoList);

	QList<ParameterTreeItem*> scalapackItems;
	scalapackItems.append(new BoolTreeItem("2D distribution", this, true));
	scalapackItems.append(new IntTreeItem("grid_nb", "grid_nb:", this, 32));
	ParameterListTreeItem *scalapackList = new ParameterListTreeItem("Scalapack", this, scalapackItems);
	addChild(scalapackList);

	QList<ParameterTreeItem*> superludistItems;
	superludistItems.append(new BoolTreeItem("Equil", this));
	superludistItems.append(new BoolTreeItem("PrintNonzeros", this));
	superludistItems.append(new BoolTreeItem("ReplaceTinyPivot", this, true));
	superludistItems.append(new BoolTreeItem("ReuseSymbolic", this));
	superludistItems.append(new IntTreeItem("perm_c", "perm_c:", this));
	superludistItems.append(new IntTreeItem("perm_r", "perm_r:", this));
	QStringList colpermSelections = QStringList() << "MMD At Plus A" << "Natural" << "MMD At A" << "COLAMD" << "My Permc";
	superludistItems.append(new StringTreeItem("ColPerm", "Column Perm:", colpermSelections, this));
	QStringList factSelections = QStringList() << "Same Pattern Same Row Perm" << "DOFACT" << "Same Pattern" << "FACTORED";
	superludistItems.append(new StringTreeItem("Fact", "Factored:", factSelections, this));
	QStringList iterrefineSelections = QStringList() << "Double" << "No" << "Extra";
	superludistItems.append(new StringTreeItem("IterRefine", "Refinement:", iterrefineSelections, this));
	QStringList rowpermSelections = QStringList() << "Large Diagonal" << "Natural" << "My Permr";
	superludistItems.append(new StringTreeItem("Row Perm", "Perm:", rowpermSelections, this));
	ParameterListTreeItem *superludistList = new ParameterListTreeItem("Superludist", this, superludistItems);
	addChild(superludistList);
}

