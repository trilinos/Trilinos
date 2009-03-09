/*
 * stratrunner.hpp
 *
 *  Created on: Jan 17, 2009
 *      Author: kurtis
 */

#ifndef STRATRUNNER_HPP_
#define STRATRUNNER_HPP_
class QTextEdit;
class Epetra_Vector;
class QString;

/**
 * Runs a Stratimikos solver on a matrix file
 */
class StratRunner{
public:
	/**
	 * Actually runs Stratimikos on a matrix file.
	 *
	 * @param xmlFileName Name of the XML file containing the specifications for the solver.
	 * @param matrixFile Name of the file containing the matrix to be solved.
	 * @param outputDisplay Where the output of the solver should go.
	 */
	static int runStrat(QString xmlFileName, QString matrixFile, QTextEdit *outputDisplay);
	/**
	 * Helper function to compute a single norm for a vector.
	 *
	 * @param Epetra_Vector vector for which to compute the norm.
	 */
	static double epetraNorm2( const Epetra_Vector &v );
};

#endif /* STRATRUNNER_HPP_ */
