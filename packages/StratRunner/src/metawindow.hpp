/*
 * metawindow.hpp
 *
 *  Created on: Nov 15, 2008
 *      Author: klnusbau
 */

#ifndef METAWINDOW_HPP_
#define METAWINDOW_HPP_
#include <QMainWindow>
class QAction;
class QMenu;
class SolverTree;

/**
 * The Main Window that contains all other widgets in the stratimikosGUI
 */
class MetaWindow : public QMainWindow{
	Q_OBJECT
public:
	/**
	 * Constructs a MainWindow object.
	 */
	MetaWindow(QString fileName=QString());
	~MetaWindow();
	friend class StartupWidget;
protected:
	void closeEvent(QCloseEvent *event);
private slots:
	/**
	 * Creates a new stratimikos solver
	 */
	void newSolve();
	/**
	 * Creates a new Run Window
	 */
	void newRunWindow();
	/**
	 * Saves the current solver settings to a user specified file.
	 */
	bool saveSolveAs();
	/**
	 * Saves the current solver to the file the user has already specified.
	 */
	void saveSolve();
	/**
	 * Loads a solver the user was previously working on and had saved.
	 */
	void loadSolve();
	/**
	 * Asks the user whether or not they would like to currently save the file they are working on. Should be used when the user has modified the file and is about to perform an action that would cause those modifiation to be lost.
	 */
	bool saveCurrentUnsavedFile();
	/**
	 * Loads a document from the set of recent documents
	 */
	void loadRecentDoc();
	/**
	 * Shows information about the program.
	 */
	void showAbout();
private:
	SolverTree *theSolverTreeWidget;
	QAction *newAct, *openRunWindowAct, *loadAct, *saveAct, *saveAsAct, *quitAct, *aboutAct;
	QMenu *fileMenu, *runMenu, *recentMenu, *helpMenu;
	QString currentLoadDir, currentSaveDir;
	QStringList recentDocsList;
	void createMenus();
	void createActions();
	void load();
	void loadLastSettings();
	void saveSettings();
	void addRecentDocument(QString recentDocument);
	void updateRecentDocsMenu();
};

#endif /* METAWINDOW_HPP_ */
