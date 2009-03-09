/*
 * runwindow.hpp
 *
 *  Created on: Jan 17, 2009
 *      Author: kurtis
 */

#ifndef RUNWINDOW_HPP_
#define RUNWINDOW_HPP_
#include <QWidget>
class QString;
class QTextEdit;
class QLineEdit;
class QPushButton;
class QGridLayout;

/**
 * The RunWindow widget is used for running a solver that has been created. It uses the current solver being used in the SolverTree and feeds in a matrix supplied by the user. The results are displayed in an output window.
 */
class RunWindow : public QWidget{
	Q_OBJECT
public:
	/**
	 * Constructs a RunWindow object.
	 *
	 * @param xmlFileName Name of the xml file associated with the current solver tree.
	 */
	RunWindow(QString xmlFileName);
private slots:
	/**
	 * Opens a matrix file to be solved.
	 */
	void openMatrixFile();
	/**
	 * Runs the solver on the specified matrix file.
	 */
	void runSolver();
	/**
	 * Saves the output of the current solver
	 */
	void saveOutput();
private:
	bool solverWasRun, currentlyRunning;
	QTextEdit *displayArea;	
	QLineEdit *matrixFilePath;
	QPushButton *openMatrixFileButton, *runButton, *saveOutputButton;
	QGridLayout *theLayout;
	QString xmlFileName;
};





#endif /* RUNWINDOW_HPP_ */
