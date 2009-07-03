#include "TivaBuena_treeview.hpp"
#include <QMessageBox>
namespace TivaBuena{


TreeView::TreeView(TreeModel *treeModel, Delegate *delegate):QTreeView(){
	setModel(treeModel);
	setItemDelegateForColumn(1, delegate);
	if(treeModel->hasDependencies()){
		connect(treeModel, SIGNAL(hideData(int, const QModelIndex&)), this, SLOT(hideRow(int, const QModelIndex&)));
		connect(treeModel, SIGNAL(showData(int, const QModelIndex&)), this, SLOT(showRow(int, const QModelIndex&)));
		connect(treeModel, SIGNAL(badValue(QModelIndex, QString)), this, SLOT(handleBadValue(QModelIndex, QString)));
		connect(delegate, SIGNAL(closeEditor(QWidget*, QAbstractItemDelegate::EndEditHint)), this, SLOT(checkForOtherBadValues()));
		treeModel->issueInitilizationSignals();
	}
}

void TreeView::showRow(int row, const QModelIndex& parent){
	if(isRowHidden(row, parent)){
		setRowHidden(row, parent, false);
	}
}

void TreeView::hideRow(int row, const QModelIndex& parent){
	if(!isRowHidden(row, parent)){
		setRowHidden(row, parent, true);
	}
}

void TreeView::handleBadValue(QModelIndex badValueIndex, QString message){
	if(state() != EditingState && !isRowHidden(badValueIndex.row(), badValueIndex.parent())){
		QMessageBox::warning(this, "Bad parameter value", message);
		setCurrentIndex(badValueIndex);
		edit(badValueIndex);
	}
	else if(!isRowHidden(badValueIndex.row(), badValueIndex.parent())){
		invalidInicies.enqueue(invalidIndex(badValueIndex, message));
	}
}

void TreeView::checkForOtherBadValues(){
	if(invalidInicies.size() != 0){
		invalidIndex needsToBeEdited = invalidInicies.dequeue();
		handleBadValue(needsToBeEdited.first, needsToBeEdited.second);
	}
}





}

