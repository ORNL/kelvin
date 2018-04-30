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
 * The mass can be retrieved as elements m_{ij} or as a
 */
class MassMatrix {
public:
	MassMatrix();
	virtual ~MassMatrix();
};

} /* namespace Kelvin */

#endif /* SRC_MASSMATRIX_H_ */
