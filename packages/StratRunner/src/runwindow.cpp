#include "runwindow.hpp"
#include "stratrunner.hpp"
#include <QTextEdit>
#include <QGridLayout>
#include <QFileDialog>
#include <QPushButton>
#include <QLineEdit>
#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include <QTextStream>

RunWindow::RunWindow(QString xmlFileName)
	:QWidget(0)
{
	setAttribute(Qt::WA_DeleteOnClose, true);
	this->xmlFileName = xmlFileName;	
	displayArea = new QTextEdit(this);
	displayArea->setReadOnly(true);
	solverWasRun = false;
	currentlyRunning = false;
	lastSaveDir = "";

	openMatrixFileButton = new QPushButton(tr("Open Matrix File"), this);
	connect(openMatrixFileButton, SIGNAL(clicked()), this, SLOT(openMatrixFile()));
	runButton = new QPushButton(tr("Run"), this);
	connect(runButton, SIGNAL(clicked()), this, SLOT(runSolver()));
	matrixFilePath = new QLineEdit(tr(""), this);
	saveOutputButton = new QPushButton(tr("Save Output"), this);
	connect(saveOutputButton, SIGNAL(clicked()), this, SLOT(saveOutput()));
	theLayout = new QGridLayout();
	theLayout->addWidget(openMatrixFileButton, 0,0);
	theLayout->addWidget(matrixFilePath,0,1);
	theLayout->addWidget(runButton,1,0);
	theLayout->addWidget(displayArea,2,0,1,3);
	theLayout->addWidget(saveOutputButton,3,2);
	setLayout(theLayout);
	resize(640,480);
	loadLastSettings();

}

RunWindow::~RunWindow(){
	saveSettings();
}

void RunWindow::saveSettings(){
	QFile *file = new QFile("runsettings.xml");
	file->open(QIODevice::WriteOnly);
	QXmlStreamWriter xmlWriter(file);

	xmlWriter.setAutoFormatting(true);
	xmlWriter.writeStartDocument();
	xmlWriter.writeStartElement("settings");
		xmlWriter.writeStartElement("lastloaddir");
		QString lastLoadPath = matrixFilePath->text();
		if(lastLoadPath != ""){
			xmlWriter.writeCharacters(lastLoadPath);
		}
		xmlWriter.writeEndElement();
		xmlWriter.writeStartElement("lastsavedir");
		xmlWriter.writeCharacters(lastSaveDir);
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
	xmlWriter.writeEndElement();
	xmlWriter.writeEndDocument();

	file->close();
	delete file;
}

void RunWindow::loadLastSettings(){
	QFile *file = new QFile("runsettings.xml");
	if(file->open(QIODevice::ReadOnly)){
		QXmlStreamReader xmlReader(file);
		while(!xmlReader.isEndDocument()){
			if(xmlReader.isStartElement()){
				if(xmlReader.name() == "lastsavedir"){
					QString tempSave = xmlReader.readElementText();
					if(tempSave != "")
						lastSaveDir = tempSave;
				}
				else if(xmlReader.name() == "lastloaddir"){
					QString tempCurLoad = xmlReader.readElementText();
					if(tempCurLoad != "")
						matrixFilePath->setText(tempCurLoad);
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
			}
			xmlReader.readNext();
		}
		file->close();
	}
	delete file;
}


		


void RunWindow::openMatrixFile(){
	matrixFilePath->setText(QFileDialog::getOpenFileName(this, tr("Open Matrix File"), QDir::homePath(), tr("Matrix Files (*.mtx)")));
}	

void RunWindow::runSolver(){
	currentlyRunning = true;
	if(!matrixFilePath->text().isEmpty()){
		StratRunner::runStrat(xmlFileName, matrixFilePath->text(), displayArea);		
	}
	else{
		QMessageBox::critical(this, tr("No Matrix File"), tr("Couldn't run the solver because a matrix file was not specified"));
	}
	solverWasRun = true;
	currentlyRunning = false;
}

void RunWindow::saveOutput(){
	if(solverWasRun && !currentlyRunning){
		QString saveFileName = QFileDialog::getSaveFileName(this, QString(tr("Save output...")), QString(lastSaveDir), tr("Text (*.txt)"));
		if(saveFileName != ""){
			QFile saveFile(saveFileName);
			if(!saveFile.open(QIODevice::WriteOnly | QIODevice::Text)){
				QMessageBox::critical(this, tr("Can't open file"), tr("Sorry, but I can't seem to open that file"));
				return;
			}

			QTextStream out(&saveFile);
			out << displayArea->toPlainText();
			
			saveFile.close();
			lastSaveDir = saveFileName;
		}
	}
	else{
		if(currentlyRunning){
			QMessageBox::critical(this, tr("Currently Running"), tr("Hey there. Hold your horses. The current solver is still running"));
		}
		else if(!solverWasRun){
			QMessageBox::critical(this, tr("No solver run"), tr("You'll need to run a solver first before you can save any output."));
		}
	}
}

