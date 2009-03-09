#include "startupwidget.hpp"
#include "metawindow.hpp"
#include <QGridLayout>
#include <QPushButton>
#include <QLabel>


StartupWidget::StartupWidget(MetaWindow *theMeta)
	:QDialog(theMeta)
{
	setWindowModality(Qt::ApplicationModal);
	setAttribute(Qt::WA_DeleteOnClose, true);
	this->theMeta = theMeta;

	instructionLabel = new QLabel(tr("Would you like to load a previously saved solver, or start a new one?"), this);
	loadFileButton =  new QPushButton(tr("Load Solver"), this);
	connect(loadFileButton, SIGNAL(clicked()), this, SLOT(loadSlot()));
	newSolverButton = new QPushButton(tr("New Solver"), this);
	connect(newSolverButton, SIGNAL(clicked()), this, SLOT(newSolverSlot()));
	theLayout = new QGridLayout();
	theLayout->addWidget(instructionLabel, 0,0,1,2);
	theLayout->addWidget(loadFileButton,1,0);
	theLayout->addWidget(newSolverButton,1,1);
	setLayout(theLayout);
}

void StartupWidget::loadSlot(){
	theMeta->load();
	close();
}

void StartupWidget::newSolverSlot(){
	close();
}

