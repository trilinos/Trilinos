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
#include <QErrorMessage>
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
		QString saveFileName = QFileDialog::getSaveFileName(this, QString(tr("Save output...")), QString(QDir::homePath() + "/SolverOutput.txt"), tr("Text (*.txt)"));
		if(saveFileName != ""){
			QFile saveFile(saveFileName);
			if(!saveFile.open(QIODevice::WriteOnly | QIODevice::Text)){
				QErrorMessage *unableToOpenError = new QErrorMessage(this);
				unableToOpenError->showMessage("Unable to open that file");
				return;
			}

			QTextStream out(&saveFile);
			out << displayArea->toPlainText();
			
			saveFile.close();
		}
	}
	else{
		QErrorMessage *noSolverRun = new QErrorMessage(this);
		if(currentlyRunning){
			noSolverRun->showMessage(QString(tr("Please wait for the solver you are currently running to finish")));
		}
		else if(!solverWasRun){
			noSolverRun->showMessage(QString(tr("You have not run a solver yet. Please run a solver first")));
		}
	}
}

