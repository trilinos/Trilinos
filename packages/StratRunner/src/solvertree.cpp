/*
 * TODO
 * -(maybe) all lists (QStringlists, list of tree items, etc.) should just be in alphabetical order
 */
#include "solvertree.hpp"
#include "stratroot.hpp"
#include <QTreeWidgetItem>
#include <QStringList>
#include <QInputDialog>
#include <QFileDialog>
#include <QMenu>

SolverTree::SolverTree(QWidget *parent, QString saveFileName)
	:QTreeWidget(parent)
{
	setContextMenuPolicy(Qt::CustomContextMenu);
	connect(this, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(displayContextMenu(const QPoint&)));
	setColumnCount(3);
	solverRoot = new StratRoot(this);
	addTopLevelItem(solverRoot);

	QStringList headers;
	headers.append(tr("Parameter"));
	headers.append(tr("Value"));
	headers.append(tr("Type"));
	setHeaderLabels(headers);
	
	this->saveFileName = saveFileName;
	if(saveFileName != ""){
		saved = true;
		readInput(saveFileName);
	}
	else{
		saved = false;
	}
	connect(this, SIGNAL(itemChanged(QTreeWidgetItem*,int)), this, SLOT(currentFileNowModified()));
	resizeColumnToContents(0);
}

void SolverTree::displayContextMenu(const QPoint & pos){
	QStringList optionsList = dynamic_cast<ParameterTreeItem*>(currentItem())->getMenuOptions();

	if(optionsList.size() != 0){
		QMenu *contextMenu = new QMenu();
		QStringList::const_iterator optionsIt;
		for(optionsIt = optionsList.constBegin(); optionsIt != optionsList.constEnd(); optionsIt++){
			QAction *newOption = new QAction(*optionsIt, this);
			contextMenu->addAction(newOption);
		}
		connect(contextMenu, SIGNAL(triggered(QAction*)), this, SLOT(launchDialog(QAction*)));
		contextMenu->popup(mapToGlobal(pos));
	}

}

void SolverTree::launchDialog(QAction *option){
	QString optionName = option->text();
	ParameterTreeItem *curItem = dynamic_cast<ParameterTreeItem*>(currentItem());
	QString optionType = curItem->getOptionType(optionName);
	QStringList optionOptionsList = curItem->getOptionOptions(optionName);
	bool ok;
	QString result;
	if(optionType == "List"){
		result = QInputDialog::getItem(this, optionName, curItem->getOptionLabel(optionName), optionOptionsList, optionOptionsList.indexOf(curItem->text(1)), false, &ok);
	}
	else if(optionType == "bool"){
		ok = true;
		if(optionName == "Set True"){
			result = "true";
		}
		else{
			result = "false";
		}
	}
	else if(optionType == "int"){
		result = QString::number(QInputDialog::getInteger(this, optionName, curItem->getOptionLabel(optionName), optionOptionsList.at(0).toInt(), optionOptionsList.at(1).toInt(), optionOptionsList.at(2).toInt(), 1, &ok)); 
	}
	else if(optionType == "double"){
		result = QString::number(QInputDialog::getDouble(this, optionName, curItem->getOptionLabel(optionName), optionOptionsList.at(0).toDouble(), optionOptionsList.at(1).toDouble(), optionOptionsList.at(2).toDouble(), 1, &ok)); 
	}
	else if(optionType == "filename"){
		ok = true;	
		result = QFileDialog::getSaveFileName(this, curItem->getOptionLabel(optionName), optionOptionsList.at(0));
	}
	else if(optionType =="freestring"){
		result = QInputDialog::getText(this, optionName, curItem->getOptionLabel(optionName), QLineEdit::Normal, curItem->text(1), &ok);
	}
	if(ok && !result.isEmpty()){
		curItem->changeParameter(optionName, result);
	}
}

bool SolverTree::writeOutput(QString fileName){
	QFile *file = new QFile(fileName);
	if(!file->open(QIODevice::WriteOnly))
		return false;
	QXmlStreamWriter xmlWriter(file);
	xmlWriter.setAutoFormatting(true);
	xmlWriter.writeStartDocument();
	solverRoot->writeOutput(xmlWriter);
	xmlWriter.writeEndDocument();
	file->close();
	delete file;
	saved = true;
	saveFileName = fileName;
	return true;
}

void SolverTree::readInput(QString fileName){
	QFile *file = new QFile(fileName);
	file->open(QIODevice::ReadOnly);
	QXmlStreamReader xmlReader(file);
	solverRoot->readInput(xmlReader);
	file->close();
	delete file;
	saved = true;
	saveFileName = fileName;
}

QString SolverTree::getSaveFileName(){
	return saveFileName;
}

bool SolverTree::isSaved(){
	return saved;
}

void SolverTree::currentFileNowModified(){
	saved = false;
}

void SolverTree::reset(){
	delete solverRoot;
	solverRoot = new StratRoot(this);
	addTopLevelItem(solverRoot);
	saveFileName = "";
	currentFileNowModified();
}

