/*
 * mlsettings.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
/*
 * TODO
 * -find list for smoother: type
 */
#include "mlsettingsitem.hpp"
#include "booltreeitem.hpp"
#include "doubletreeitem.hpp"
#include "inttreeitem.hpp"
#include "stringtreeitem.hpp"

MLSettingsItem::MLSettingsItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, MLSettingsItem::Type)
{
	setText(0,QString("ML Settings"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void MLSettingsItem::addParameters(){
	addChild(new DoubleTreeItem("aggregation: damping factor", "Factor:", this, 1.333));
	addChild(new DoubleTreeItem("aggregation: edge prolongator drop threshold", "Threshold:", this));
	addChild(new IntTreeItem("aggregation: local aggregates", "value:", this, 1));
	addChild(new IntTreeItem("aggregation: next-level aggregates per process", "value:", this, 128));
	addChild(new IntTreeItem("aggregation: nodes per aggregate", "value:", this, 512));
	addChild(new IntTreeItem("coarse: max size", "Size:", this, 128));
	addChild(new IntTreeItem("eigen-analysis: iterations", "Iteration:", this, 10));
	addChild(new IntTreeItem("max levels", "Max:", this, 10));
	addChild(new BoolTreeItem("smoother: Aztec as solver", this));
	addChild(new BoolTreeItem("smoother: Hiptmair efficient symmetric", this, true));
	addChild(new DoubleTreeItem("smoother: damping factor", "Factor:", this, 1));
	
	QStringList smootherPrePostList = QStringList() << "Both" << "Pre" << "Post";
	addChild(new StringTreeItem("smoother: pre or post", "Pre Or Post:", smootherPrePostList, this));
	QStringList subSmootherTypeList = QStringList() << "Chebyshev"<< "symmetric Gauss-Seidel"<< "ILU"<< "IC"<< "ILUT"<< "ICT";
	addChild(new StringTreeItem("subsmoother: type", "Type:", subSmootherTypeList,this));
	QStringList aggregationTypeList = QStringList() << "Uncoupled-MIS"<< "MIS"<< "METIS"<< "Uncoupled"<< "ParMETIS";
	addChild(new StringTreeItem("aggregation: type", "Type:",  aggregationTypeList, this));
	QStringList coarseTypeList = QStringList() << "Amesos-KLU"<< "Jacobi"<< "Gauss-Seidel"<< "ML Gauss-Seidel"<< "symmetric Gauss-Seidel"<< "ML symmetric Gauss-Seidel"<< "Chebyshev"<< "Hiptmair"<< "Amesos-Superlu"<< "Amesos-UMFPACK"<< "Amesos-Superludist"<< "Amesos-MUMPS"<< "user-defined"<< "SuperLU";
	addChild(new StringTreeItem("coarse: type", "Type:", coarseTypeList, this));
	QStringList defaultValuesList = QStringList() << "maxwell" << "SA" << "DD" << "DD-ML" << "NSSA";
	addChild(new StringTreeItem("default values", "Default Values:", defaultValuesList, this));
	QStringList eigenAnalysisTypeList = QStringList() << "cg" << "Anorm" << "power-method";
	addChild(new StringTreeItem("eigen-analysis: type", "Type:", eigenAnalysisTypeList, this));
	QStringList incDecList = QStringList() << "increasing" << "decreasing";
	addChild(new StringTreeItem("increasing or decreasing", "Increasing or Decreasing:", incDecList, this));
	QStringList precTypeList = QStringList() << "MGV" << "MGW" << "full-MGV" << "one-level-postsmoothing" << "two-level-additive" << "two-level-hybrid" << "two-level-hybrid2" << "projected MGV";
	addChild(new StringTreeItem("prec type", "Type:", precTypeList, this));
	QStringList smootherTypeList = QStringList() << "Chebyshev" << "symmetric Gauss-Seidel";
	addChild(new StringTreeItem("smoother: type", "Type:", smootherTypeList, this));

	addChild(new IntTreeItem("smoother: sweeps", "Sweeps:", this, 1));
	addChild(new DoubleTreeItem("subsmoother: Chebyshev alpha", "Value:", this, 20));
	addChild(new DoubleTreeItem("subsmoother: edge sweeps", "Value:", this, 4));
	addChild(new DoubleTreeItem("subsmoother: node sweeps", "Value:", this, 4));
}

