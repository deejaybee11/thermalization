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
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mkl.h"
#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/SaveData.hpp"


Potential::Potential(SimulationData &sim_data) {
	//Allocate memory for the arrays

	this->harmonic_trap = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->non_linear = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->kinetic_energy = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);
	this->pos_time_evolution = (MKL_Complex16*)mkl_malloc(sim_data.get_N() * sizeof(MKL_Complex16), 64);
	this->mom_time_evolution = (MKL_Complex16*)mkl_malloc(sim_data.get_N() * sizeof(MKL_Complex16), 64);

	
	//Fill arrays with values
	double h_pot_val = 0;
	double k_en_val = 0;
	int index;

	#pragma omp parallel for private(index, h_pot_val)
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			for (int k = 0; k < sim_data.get_num_z(); ++k) {

				index = i * sim_data.get_num_y() * sim_data.get_num_z() + j * sim_data.get_num_z() + k;
				h_pot_val = 0.5 * (pow(sim_data.gamma_x, 2) * pow(sim_data.x[i], 2) + pow(sim_data.gamma_y, 2) * pow(sim_data.y[j], 2) + pow(sim_data.gamma_z, 2) * pow(sim_data.z[k], 2));
				k_en_val = 0.5 * (pow(sim_data.px[i], 2.0) + pow(sim_data.py[j], 2.0) + pow(sim_data.pz[k], 2.0));

				this->harmonic_trap[index] = h_pot_val; 
				this->kinetic_energy[index] = k_en_val;
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

void Potential::assign_position_time_evolution(SimulationData &sim_data, WaveFunction &psi, bool trap_on, bool is_real) {
	
	double theta;
	if (is_real && trap_on) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = (this->non_linear[i] + this->harmonic_trap[i]) * 0.5 * sim_data.get_dt();	
			this->pos_time_evolution[i].real = cos(theta);
			this->pos_time_evolution[i].imag = -1.0 * sin(theta);	
		}
	}
	else if (is_real && !trap_on) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = this->non_linear[i]  * 0.5 * sim_data.get_dt();	
			this->pos_time_evolution[i].real = cos(theta);
			this->pos_time_evolution[i].imag = -1.0 * sin(theta);	
		}

	}
	else {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = (this->non_linear[i] + this->harmonic_trap[i]) * 0.5 * sim_data.get_dt();	
			this->pos_time_evolution[i].real = exp(-1.0 * theta);
			this->pos_time_evolution[i].imag = 0;	
		}
	}
}

void Potential::assign_momentum_time_evolution(SimulationData &sim_data, WaveFunction &psi, bool is_real) {
	
	double theta;
	if (is_real) {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = this->kinetic_energy[i] * sim_data.get_dt();	
			this->mom_time_evolution[i].real = cos(theta);
			this->mom_time_evolution[i].imag = -1.0 * sin(theta);	
		}
	}
	else {
		#pragma omp parallel for private(theta)
		for (int i = 0; i < sim_data.get_N(); ++i) {
			theta = this->kinetic_energy[i] * sim_data.get_dt();	
			this->mom_time_evolution[i].real = exp(-1.0 * theta);
			this->mom_time_evolution[i].imag = 0;	
		}
	}
}



Potential::~Potential() { 

	mkl_free(harmonic_trap);
	mkl_free(non_linear);
	mkl_free(kinetic_energy);
	mkl_free(pos_time_evolution);
	mkl_free(mom_time_evolution);

}

