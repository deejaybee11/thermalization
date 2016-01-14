/**The MIT License (MIT)
*
*Copyright (c) 2016 Dylan
*
*Permission is hereby granted, free of charge, to any person obtaining a copy
*of this software and associated documentation files (the "Software"), to deal
*in the Software without restriction, including without limitation the rights
*to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the Software is
*furnished to do so, subject to the following conditions:
*
*The above copyright notice and this permission notice shall be included in all
*copies or substantial portions of the Software.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*SOFTWARE.
*/

#ifndef THERMALIZATION_SIM_DATA_H
#define THERMALIZATION_SIM_DATA_H

#include <stdlib.h>

#include "mkl.h"

/**
 * SimulationData class
 *
 * Stores data and arrays used in simulation
 */

class SimulationData {
public:
	SimulationData(int num_x, int num_y, int num_z);
	~SimulationData();

	double *x;
	double *y;
	double *z;
	double *px;
	double *py;
	double *pz;

	double length_x;
	double length_y;
	double length_z;

	double dx;
	double dy;
	double dz;
	double dpx;
	double dpy;
	double dpz;

	double sigma_x;
	double sigma_y;
	double sigma_z;

	double gamma_x;
	double gamma_y;
	double gamma_z;

	double beta;

	double *chemical_potential;
	
	int num_R_steps;
	int num_I_steps;
	const char *R = "REAL";
	const char *I = "IMAG";

	int get_num_x() { return this->num_x; };
	int get_num_y() { return this->num_y; };
	int get_num_z() { return this->num_z; };
	int get_N() { return this->N; };

private:
	int num_x;
	int num_y;
	int num_z;
	int N;

	double dt;

};

#endif    //    THERMALIZATION_SIM_DATA_H
