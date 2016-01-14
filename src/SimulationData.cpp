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

#include "../include/SimulationData.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

//Class constructor
SimulationData::SimulationData(int num_x, int num_y, int num_z) {

	this->num_x = num_x;
	this->num_y = num_y;
	this->num_z = num_z;
	this->N = num_x * num_y * num_z;

	this->num_I_steps = 100000;
	this->num_R_steps = 100000;

	this->sigma_x = 1;
	this->sigma_y = 1;
	this->sigma_z = 1;

	this->gamma_x = 1;
	this->gamma_y = 1.2;
	this->gamma_z = 3;

	this->beta = 1;

	//Allocate memory and assign numbers to the grids
	this->x = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->y = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
	this->z = (double*)mkl_malloc(this->num_z * sizeof(double), 64);

	this->px = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
	this->py = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
	this->pz = (double*)mkl_malloc(this->num_z * sizeof(double), 64);

	this->chemical_potential = (double*)mkl_malloc(this->N * sizeof(double), 64);

	for (int i = 0; i < this->num_x; ++i) {
		this->x[i] = -0.5 * this->length_x + i * this->length_x / ((double)this->num_x);
	}
	for (int i = 0; i < this->num_y; ++i) {
		this->y[i] = -0.5 * this->length_y + i * this->length_y / ((double)this->num_y);
	}
	for (int i = 0; i < this->num_z; ++i) {
		this->z[i] = -0.5 * this->length_z + i * this->length_z / ((double)this->num_z);
	}

	this->dx = this->x[1] - this->x[0];
	this->dy = this->y[1] - this->y[0];
	this->dz = this->z[1] - this->z[0];

	//Increments for the momenum arrays
	double ax = -0.5 * this->num_x;
	double bx = 0.5 * this->num_x - 1.0;
	double ay = -0.5 * this->num_y;
	double by = 0.5 * this->num_y - 1.0;
	double az = -0.5 * this->num_z;
	double bz = 0.5 * this->num_z - 1.0;

	double step_x = (2 * M_PI / this->length_x) * ((bx - ax) / (this->num_x - 1.0));
	double step_y = (2 * M_PI / this->length_y) * ((by - ay) / (this->num_y - 1.0));
	double step_z = (2 * M_PI / this->length_z) * ((bz - az) / (this->num_z - 1.0));

	//Fill momentum arrays
	for (int i = 0; i < this->num_x; ++i) {
		this->px[i] = (2 * M_PI / this->length_x) * ax + i * step_x;
	}
	for (int i = 0; i < this->num_y; ++i) {
		this->py[i] = (2 * M_PI / this->length_y) * ay + i * step_y;
	}
	for (int i = 0; i < this->num_z; ++i) {
		this->pz[i] = (2 * M_PI / this->length_z) * az + i * step_z;
	}

};

//Class destructor
SimulationData::~SimulationData() {
	mkl_free(x);
	mkl_free(y);
	mkl_free(z);
	mkl_free(px);
	mkl_free(py);
	mkl_free(pz);
	mkl_free(chemical_potential);
};
