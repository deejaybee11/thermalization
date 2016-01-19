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

#include "../include/Solve.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdexcept>

#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/Potential.hpp"

void solve_imag(SimulationData &sim_data, WaveFunction &psi, Potential &pot_data) {
	
	double temp_real = 0;
	double temp_imag = 0;
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;

	MKL_LONG N[3]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y(); N[2] = sim_data.get_num_z();	
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, N);
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, (1.0 / (N[0] * N[1] * N[2])));
	status = DftiCommitDescriptor(handle);


	for (int i = 0; i < sim_data.num_I_steps; ++i) {

		if (i % 500 == 0) {
			std::cout << "Imaginary step " << i << " out of " << sim_data.num_I_steps << "." << std::endl;
		}

		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);

		status = DftiComputeForward(handle, psi.psi, psi.psi);

		vzMul(sim_data.get_N(), psi.psi, pot_data.mom_time_evolution, psi.psi);

		status = DftiComputeBackward(handle, psi.psi, psi.psi);

		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);

		//psi.get_abs(sim_data.get_N());
		//psi.get_norm(sim_data);

		
	}
	



}

void solve_real(SimulationData &sim_data, WaveFunction &psi, Potential &pot_data) {
	
	double temp_real = 0;
	double temp_imag = 0;
	MKL_LONG status = 0;
	DFTI_DESCRIPTOR_HANDLE handle = 0;

	MKL_LONG N[3]; N[0] = sim_data.get_num_x(); N[1] = sim_data.get_num_y(); N[2] = sim_data.get_num_z();	
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, N);
	status = DftiSetValue(handle, DFTI_BACKWARD_SCALE, (1.0 / (N[0] * N[1] * N[2])));
	status = DftiCommitDescriptor(handle);


	for (int i = 0; i < sim_data.num_R_steps; ++i) {

		if (i % 500 == 0) {
			std::cout << "Imaginary step " << i << " out of " << sim_data.num_R_steps << "." << std::endl;
		}

		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);

		status = DftiComputeForward(handle, psi.psi, psi.psi);

		vzMul(sim_data.get_N(), psi.psi, pot_data.mom_time_evolution, psi.psi);

		status = DftiComputeBackward(handle, psi.psi, psi.psi);

		vzMul(sim_data.get_N(), psi.psi, pot_data.pos_time_evolution, psi.psi);

		//psi.get_abs(sim_data.get_N());
		//psi.get_norm(sim_data);

		
	}
	



}