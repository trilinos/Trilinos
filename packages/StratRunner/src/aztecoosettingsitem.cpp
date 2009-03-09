/*
 * aztecoosettingsitem.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
#include "aztecoosettingsitem.hpp"
#include "stringtreeitem.hpp"
#include "doubletreeitem.hpp"

AztecOOSettingsItem::AztecOOSettingsItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, AztecOOSettingsItem::Type)
{
	setText(0,QString("AztecOO Settings"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
}

void AztecOOSettingsItem::addParameters(){
	QStringList aztecPrecons = QStringList() << "none" << "ilu" << "ilut" << "Jacobi" << "Symmetric Gauss-Seidel" << "Polynomial" << "Least-squares Polynomial";
	addChild(new StringTreeItem("Aztec Preconditioner", "Preconditioner:", aztecPrecons, this));

	QStringList solvers = QStringList() << "CG" << "GMRES" << "CGS" << "TFQMR" << "BiCGStab" << "LU";
	addChild(new StringTreeItem("Aztec Solver", "Solver:", solvers, this));

	QStringList tests = QStringList() << "r0" << "rhs" << "anorm" << "no scaling" << "sol";
	addChild(new StringTreeItem("Convergence Test", "Test:", tests, this));

	QStringList orthos = QStringList() << "Classical" << "Modified";
	addChild(new StringTreeItem("Orthogonalization", "Orthogonalization:", orthos, this));

	QStringList rcmOps = QStringList() << "Enabled" << "Disabled";
	addChild(new StringTreeItem("RCM Reordering", "RCM Reordering:", rcmOps, this));

	addChild(new DoubleTreeItem("Drop Tolerance", "Tolerance:", this));
	addChild(new DoubleTreeItem("Fill Factor", "Factor:", this));
	addChild(new DoubleTreeItem("Graph Fill", "Fill:", this));
	addChild(new DoubleTreeItem("Ill-Conditioning Threshold", "Threshold:", this));
	addChild(new DoubleTreeItem("Output Frequency", "Frequence:", this));
	addChild(new DoubleTreeItem("Overlap", "Overlap:", this));
	addChild(new DoubleTreeItem("Polynomial Order", "Order:", this, 3));
	addChild(new DoubleTreeItem("Size of Krylov Subspace", "Size:", this));
	addChild(new DoubleTreeItem("Steps", "Steps:", this, 3));
}

