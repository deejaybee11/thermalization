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

#include "../include/WaveFunction.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "../include/SimulationData.hpp"
#include "../include/Potential.hpp"


WaveFunction::WaveFunction(SimulationData &sim_data, double *harmonic_trap) {

	this->psi = (MKL_Complex16*)malloc(sim_data.get_N() * sizeof(MKL_Complex16));	
	this->psi_conj = (MKL_Complex16*)malloc(sim_data.get_N() * sizeof(MKL_Complex16));	
	this->abs_psi = (double*)mkl_malloc(sim_data.get_N() * sizeof(double), 64);

	int index;

	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			for (int k = 0; k < sim_data.get_num_z(); ++k) {

				index = i * sim_data.get_num_y() * sim_data.get_num_z() + j * sim_data.get_num_z() + k;
				//this->psi[index].real = exp(-0.3 * (pow(sim_data.x[i], 2.0) + pow(sim_data.y[j], 2.0) + pow(sim_data.z[k], 2.0)));
				//this->psi[index].imag = 0;
				if (sqrt(pow(sim_data.x[i],2.0) + pow(sim_data.y[j], 2.0) + pow(sim_data.z[k], 2.0)) <= 5) {
					this->psi[index].real = 1;
					this->psi[index].imag = 0;
				}
				else {
					this->psi[index].real = 0;
					this->psi[index].imag = 0;
				}
			}
		}
	}
}

void WaveFunction::create_superposition(SimulationData &sim_data) {

	MKL_Complex16 *temp_psi;
	temp_psi = (MKL_Complex16*)mkl_malloc(sim_data.get_N() * sizeof(MKL_Complex16), 64);
	int index;
	
	#pragma omp parallel for private(index)
	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			for (int k = 0; k < sim_data.get_num_z(); ++k) {
				index = i * sim_data.get_num_y() * sim_data.get_num_z() + j* sim_data.get_num_z() + k;
				temp_psi[index].real = this->psi[index].real;
				temp_psi[index].imag = this->psi[index].imag;
				this->psi[index].real += temp_psi[index].real * cos(2 * sim_data.laser_kick * sim_data.y[j]) - temp_psi[index].imag * sin(2 * sim_data.laser_kick * sim_data.y[j]);
				this->psi[index].imag += temp_psi[index].imag * cos(2 * sim_data.laser_kick * sim_data.y[j]) + temp_psi[index].real * sin(2 * sim_data.laser_kick * sim_data.y[j]);
			}
		}
	}
}

void WaveFunction::get_norm(SimulationData &sim_data) {

	double psi_sum = 0;
	double temp_real = 0;
	double temp_imag = 0;
	double norm_fac;
	
	for (int i = 0; i < sim_data.get_N(); ++i) {
		psi_sum += this->abs_psi[i];
	}

	norm_fac = sqrt(1.0 / (psi_sum * sim_data.dx * sim_data.dy * sim_data.dz));
	
	for (int i = 0; i < sim_data.get_N(); ++i) {
		temp_real = this->psi[i].real * norm_fac;
		temp_imag = this->psi[i].imag * norm_fac;
		this->psi[i].real = temp_real;
		this->psi[i].imag = temp_imag;
	}


	
}

void WaveFunction::get_abs(int N) {
	vzAbs(N, this->psi, this->abs_psi);
	vdMul(N, this->abs_psi, this->abs_psi, this->abs_psi);
}

WaveFunction::~WaveFunction() {
	
	mkl_free(psi);
	mkl_free(abs_psi);

}
