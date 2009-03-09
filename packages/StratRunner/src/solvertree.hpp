/*
 * solvertree.hpp
 *
 *  Created on: Dec 3, 2008
 *      Author: kurtis
 */

#ifndef SOLVERTREE_HPP_
#define SOLVERTREE_HPP_
#include <QTreeWidget>
class QTreeWidgetItem;
class StratRoot;

/**
 * A QTreeWidget that contains only ParameterTreeItems. It can display context menus from any of the ParameterTreeItems.
 */
class SolverTree : public QTreeWidget{
	Q_OBJECT
public:
	/**
	 * Constructs a SolverTree object.
	 *
	 * @param parent The parent widget.
	 */
	SolverTree(QWidget *parent=0, QString saveFileName="");
	/**
	 * Writes the xml file to be fed into stratimkios based on the input specified by the user.
	 *
	 * @param fileName The name of the file to which the SolverTree should write the XML output.
	 */
	void writeOutput(QString fileName);
	/**
	 * Reads an xml file that specifies a stratimikos solver and sets all the parameters accordingly
	 *
	 * @param fileName The name of the file from which the SolverTree should read parameter values.
	 */
	void readInput(QString fileName);
	/**
	 * Gets the name of the save file with which this SolverTree widget is associated. If this SolverTree widget has yet to be saved and thus has no save file associated with it, the funtion will return an empty string.
	 *
	 * @return The name of the save file with which this SolverTree widget is associated.
	 */
	QString getSaveFileName();
	/**
	 * Determines wether or not the current state of SolverTree has been saved.
	 *
	 * @return True if the current state of the SolverTree has been saved. False otherwise.
	 */
	bool isSaved();
	/**
	 * Resets the state of the SolverTree and all of its children to their defaults.
	 */
	void reset();
private slots:
	/**
	 * Displays the context menu for the selected ParameterTreeItem.
	 *
	 * @param pos Where the context menu should be displayed.
	 */
	void displayContextMenu(const QPoint & pos);
	/**
	 * Lauches a QDialog. The type of QDialog launched will depend on the QAction specified by the option parameter.
	 *
	 * @param option What type of QDialog should be launched.
	 */
	void launchDialog(QAction *option);
	/**
	 * When the state of any of the SolverTree's items is changed, this slot should be called
	 */
	void currentFileNowModified();
private:
	QString saveFileName;
	bool saved;
	StratRoot *solverRoot;
};

#endif /* SOLVERTREE_HPP_ */
