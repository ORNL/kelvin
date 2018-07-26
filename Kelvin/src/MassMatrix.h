/*
 * MassMatrix.h
 *
 *  Created on: Apr 30, 2018
 *      Author: bkj
 */

#ifndef SRC_MASSMATRIX_H_
#define SRC_MASSMATRIX_H_

namespace Kelvin {

/**
 * The mass matrix represents the mass of the system discretized on the
 * computational grid. A mass matrix element m_{ij} represents the mass shared
 * between nodes i and j.
 *
 * The mass can be retrieved as matrix elements m_{ij} or as a list of diagonal
 * entries on mass-lumped diagonalized matrix, m_{D,ij}. For descriptions of
 * both see Sulsky's 1994 paper "A particle method for history-dependent
 * materials."
 */
class MassMatrix {
public:
	MassMatrix();
	virtual ~MassMatrix();
};

} /* namespace Kelvin */

#endif /* SRC_MASSMATRIX_H_ */
