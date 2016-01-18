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

#include "../include/Potential.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mkl.h"
#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"

Potential::Potential(SimulationData &sim_data) {
	//Allocate memory for the arrays
	
	this->harmonic_trap = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->non_linear = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->time_evolution = (MKL_Complex16*)mkl_malloc(sim_data.get_N() * sizeof(MKL_Complex16), 64);
	
	//Fill arrays with values
	double h_pot_val = 0;
	int index;

	#pragma omp parallel for private(index, h_pot_val)
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			for (int k = 0; k < sim_data.get_num_z(); ++k) {

				index = i * sim_data.get_num_x() * sim_data.get_num_z() + j * sim_data.get_num_z() + k;
				h_pot_val = 0.5 * (pow(sim_data.gamma_x, 2) * pow(sim_data.x[i], 2) + pow(sim_data.gamma_y, 2) * pow(sim_data.y[j], 2) + pow(sim_data.gamma_z, 2) * pow(sim_data.z[k], 2));

				this->harmonic_trap[index] = h_pot_val; 

			}
		}
	}
}

void Potential::calculate_non_linear_energy(SimulationData &sim_data, WaveFunction &psi) {

	double non_linear_val;

	#pragma omp parallel for private(non_linear_val)
	for (int i = 0; i < sim_data.get_N(); ++i) {
		non_linear_val = sim_data.beta * psi.abs_psi[i];
		this->non_linear[i] = non_linear_val;
	}
}

void Potential::assign_time_evolution(SimulationData &sim_data, WaveFunction &psi, Potential &potential_data, bool trap_on, bool is_real) {
	
	double theta;
	#pragma omp parallel for private(theta)
	for (int i = 0; i < sim_data.get_N(); ++i) {
		theta = (this->non_linear[i] + this->harmonic_trap[i]) * sim_data.get_dt();	
		time_evolution[i].real = cos(theta);
		time_evolution[i].imag = -1.0 * sin(theta);	
	}
}

Potential::~Potential() {

	mkl_free(harmonic_trap);
	mkl_free(non_linear);
}

