/*
 * startupwidget.hpp
 *
 *  Created on: Feb 15, 2009
 *      Author: kurtis
 */

#ifndef STARTUPWIDGET_HPP_
#define STARTUPWIDGET_HPP_
#include <QDialog>
class QPushButton;
class QGridLayout;
class MetaWindow;
class QLabel;

/**
 * The StartupWidget widget give the user and initial choice as to whether they would like to start a new solver or load an old one.
 */
class StartupWidget : public QDialog{
	Q_OBJECT
public:
	/**
	 * Constructs a StartupWidget object.
	 *
	 */
	StartupWidget(MetaWindow *theMeta);
private slots:
	void loadSlot();
	void newSolverSlot();
private:
	QLabel *instructionLabel;
	QPushButton *loadFileButton, *newSolverButton;
	QGridLayout *theLayout;
	MetaWindow *theMeta;
};

#endif /* STARTUPWIDGET_HPP_ */

