/*
 * metawindow.cpp
 *
 *  Created on: Nov 15, 2008
 *      Author: klnusbau
 */

#include "runwindow.hpp"
#include "metawindow.hpp"
#include "solvertree.hpp"
#include <QFileDialog>
#include <QMessageBox> 
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QtGui>
#include <QIcon>
const int numRecentDocuments = 7;

MetaWindow::MetaWindow(QString fileName){
	theSolverTreeWidget = new SolverTree(this, fileName);
	setCentralWidget(theSolverTreeWidget);
	createActions();
	createMenus();
	resize(800,600);
	currentLoadDir = QDir::homePath();
	currentSaveDir = QDir::homePath();
	loadLastSettings();
	setWindowTitle(tr("StratRunner"));
	setWindowIcon(QIcon("stratrunner.png"));

}

MetaWindow::~MetaWindow(){
	saveSettings();
}

void MetaWindow::saveSettings(){
	QFile *file = new QFile("settings.xml");
	file->open(QIODevice::WriteOnly);
	QXmlStreamWriter xmlWriter(file);

	xmlWriter.setAutoFormatting(true);
	xmlWriter.writeStartDocument();
	xmlWriter.writeStartElement("settings");
		xmlWriter.writeStartElement("lastsavedir");
		xmlWriter.writeCharacters(currentSaveDir);
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("lastloaddir");
		xmlWriter.writeCharacters(currentLoadDir);
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("xres");
		xmlWriter.writeCharacters(QString::number(width()));
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("yres");
		xmlWriter.writeCharacters(QString::number(height()));
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("xpos");
		xmlWriter.writeCharacters(QString::number(x()));
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("ypos");
		xmlWriter.writeCharacters(QString::number(y()));
		xmlWriter.writeEndElement();
		for(int i =0; i<recentDocsList.size(); i++){
			xmlWriter.writeStartElement("recentdoc");
				xmlWriter.writeCharacters(recentDocsList.at(i));
			xmlWriter.writeEndElement();
		}
	xmlWriter.writeEndElement();
	xmlWriter.writeEndDocument();

	file->close();
	delete file;
}
	
void MetaWindow::loadLastSettings(){
	QFile *file = new QFile("settings.xml");
	if(file->open(QIODevice::ReadOnly)){
		QXmlStreamReader xmlReader(file);
		while(!xmlReader.isEndDocument()){
			if(xmlReader.isStartElement()){
				if(xmlReader.name() == "lastsavedir"){
					QString tempCurSave = xmlReader.readElementText();
					if(tempCurSave != "")
						currentSaveDir = tempCurSave;
				}
				else if(xmlReader.name() == "lastloaddir"){
					QString tempCurLoad = xmlReader.readElementText();
					if(tempCurLoad != "")
						currentLoadDir = tempCurLoad;
				}
				else if(xmlReader.name() == "xres"){
					QString theWidth = xmlReader.readElementText();
					if(theWidth != "")
						resize(theWidth.toInt(), height());
				}
				else if(xmlReader.name() == "yres"){
					QString theHeight = xmlReader.readElementText();
					if(theHeight != "")
						resize(width(), theHeight.toInt());
				}
				else if(xmlReader.name() == "xpos"){
					QString xpos = xmlReader.readElementText();
					if(xpos != "")
						move(xpos.toInt(), y());
				}
				else if(xmlReader.name() == "ypos"){
					QString ypos = xmlReader.readElementText();
					if(ypos != "")
						move(x(), ypos.toInt());
				}
				else if(xmlReader.name() == "recentdoc"){
					addRecentDocument(xmlReader.readElementText());
				}
			}
			xmlReader.readNext();
		}
		file->close();
	}
	delete file;


}

void MetaWindow::addRecentDocument(QString recentDocument){
	recentDocsList.prepend(recentDocument);
	if(recentDocsList.size() > numRecentDocuments){
		recentDocsList.removeLast();
	}
//	updateRecentDocsMenu();
}

void MetaWindow::updateRecentDocsMenu(){
	recentMenu->clear();
	for(int i=0; i<recentDocsList.size(); i++){
		QAction *recentDocAct = new QAction(recentDocsList.at(i).section("/",-1,-1),this);
		connect(recentDocAct, SIGNAL(triggered()), this, SLOT(loadRecentDoc()));
		recentMenu->addAction(recentDocAct);
	}
}
	


void MetaWindow::createActions(){
	newAct = new QAction(tr("&New"),this);
	newAct->setShortcut(tr("Ctrl+N"));
	newAct->setStatusTip(tr("New Simulation"));
	connect(newAct, SIGNAL(triggered()), this, SLOT(newSolve()));

	openRunWindowAct = new QAction(tr("Open &Run Window"),this);
	openRunWindowAct->setShortcut(tr("Ctrl+R"));
	openRunWindowAct->setStatusTip(tr("Open a Window to Run this Stratimikos Simulation"));
	connect(openRunWindowAct, SIGNAL(triggered()), this, SLOT(newRunWindow()));

	saveAct = new QAction(tr("&Save"),this);
	saveAct->setShortcut(tr("Ctrl+S"));
	saveAct->setStatusTip(tr("Save the current file."));
	connect(saveAct, SIGNAL(triggered()), this, SLOT(saveSolve()));

	saveAsAct = new QAction(tr("Save As..."),this);
	saveAsAct->setStatusTip(tr("Save the current file to a specified file name."));
	connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveSolveAs()));

	loadAct = new QAction(tr("&Load"),this);
	loadAct->setShortcut(tr("Ctrl+L"));
	loadAct->setStatusTip(tr("Load Simulation"));
	connect(loadAct, SIGNAL(triggered()), this, SLOT(loadSolve()));

	quitAct = new QAction(tr("&Quit"),this);
	quitAct->setShortcut(tr("Ctrl+Q"));
	quitAct->setStatusTip(tr("Quit Easy Tramonto"));
	connect(quitAct, SIGNAL(triggered()), qApp, SLOT(quit()));

	aboutAct = new QAction(tr("About"),this);
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(showAbout()));

}

void MetaWindow::createMenus(){
//	recentMenu = new QMenu(tr("Recent Solvers"));
//	QAction *noRecentAct = new QAction(tr("No Recent Documents"),this);
//	noRecentAct->setEnabled(false);
//	recentMenu->addAction(noRecentAct);
	fileMenu = menuBar()->addMenu(tr("File"));
	fileMenu->addAction(newAct);
	//fileMenu->addMenu(recentMenu);
	fileMenu->addSeparator();
	fileMenu->addAction(saveAct);
	fileMenu->addAction(saveAsAct);
	fileMenu->addAction(loadAct);
	fileMenu->addSeparator();
	fileMenu->addAction(quitAct);
	runMenu = menuBar()->addMenu(tr("Run"));
	runMenu->addAction(openRunWindowAct);
	helpMenu = menuBar()->addMenu(tr("Help"));
	helpMenu->addAction(aboutAct);
}

void MetaWindow::newSolve(){
	if(!theSolverTreeWidget->isSaved()){
		saveCurrentUnsavedFile();
	}
	theSolverTreeWidget->reset();
}

void MetaWindow::newRunWindow(){
	if(!theSolverTreeWidget->isSaved()){
		QMessageBox::warning(this, QString(tr("Unsaved Changes!")), QString(tr("You have made changes to the solver that you have not yet saved. Please save them before continuing")), QMessageBox::Save, QMessageBox::Save);
		saveSolveAs();
	}
	if(theSolverTreeWidget->isSaved()){
		RunWindow *theRunWindow = new RunWindow(theSolverTreeWidget->getSaveFileName());
		theRunWindow->show();
	}
}

void MetaWindow::saveSolve(){
	QString currentFileName = theSolverTreeWidget->getSaveFileName();
	if(currentFileName != ""){
		theSolverTreeWidget->writeOutput(currentFileName);
	}
	else{
		saveSolveAs();
	}
}

bool MetaWindow::saveSolveAs(){
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save To..."), currentSaveDir);
	if(fileName != ""){
		if(theSolverTreeWidget->writeOutput(fileName)){
			currentSaveDir = fileName.section("/",0,-2);
			addRecentDocument(fileName);		
			return true;
		}
	}
	return false;
}

bool MetaWindow::saveCurrentUnsavedFile(){
		QMessageBox saveQuestion(QMessageBox::Question, tr("Save?"), tr("The current solver you are working on has not be saved. Would you like to save it now?"), QMessageBox::Yes, this);
		saveQuestion.addButton(QMessageBox::No);
		int shouldSave = saveQuestion.exec(); 
		if(shouldSave == QMessageBox::Yes){
			return saveSolveAs();
		}
		return true;
}

void MetaWindow::loadSolve(){
	if(!theSolverTreeWidget->isSaved()){
		saveCurrentUnsavedFile();
	}
	load();
}

void MetaWindow::loadRecentDoc(){
	QString docName = dynamic_cast<QAction*>(sender())->text();
	int i =0;
	for(; i<recentDocsList.size();i++){
		if(recentDocsList.at(i).contains(docName)){
			break;
		}
	}
	if(!theSolverTreeWidget->isSaved()){
		if(saveCurrentUnsavedFile()){
			theSolverTreeWidget->readInput(recentDocsList.at(i));
		}
	}
}

void MetaWindow::load(){
	QString fileName = QFileDialog::getOpenFileName(this, tr("Load..."), currentLoadDir, tr("Xml (*.xml)"));
	if(fileName != ""){
		theSolverTreeWidget->readInput(fileName);
		currentLoadDir = fileName.section("/",0,-2);
		addRecentDocument(fileName);
	}
}

void MetaWindow::showAbout(){
	QMessageBox::about(this, "About StratRunner", "StratRunner was developed by Kurtis Nusbaum. For detailed information on how to use StratRunner please consult the UsersGuide.pdf.\nLicense:LGPL.\nContact: klnusbaum@gmail.com");
}

void MetaWindow::closeEvent(QCloseEvent *event){
	if(!theSolverTreeWidget->isSaved()){
		if(saveCurrentUnsavedFile()){
			event->accept();
			saveSettings();
		}
		else{
			event->ignore();
		}
	}
	else{
		event->accept();
		saveSettings();
	}
}

