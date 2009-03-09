#include <QApplication>
#include <QtGui>
#include "metawindow.hpp"
#include "startupwidget.hpp"

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	MetaWindow *theWindow = new MetaWindow(argv[1]);
	theWindow->show();
	if(argc < 2){
		StartupWidget *theStartup = new StartupWidget(theWindow);
		theStartup->show();
	}
	return a.exec();
}
