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

#include <stdlib.h>
#include <iostream>


#include "mkl.h"

#include "../include/SimulationData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/Potential.hpp"
#include "../include/SaveData.hpp"
#include "../include/Solve.hpp"

#if !defined(MKL_ILP64)
	#define LI "%li"
	#else
	#define LI "%lli"
#endif


int main() {

	putenv("KMP_BLOCKTIME=0");
	putenv("KMP_AFFINITY=verbose,granularity=fine,compact,norespect");
	system("echo KMP_BLOCKTIME = $KMP_BLOCKTIME");
	system("echo KMP_AFFINITY = $KMP_AFFINITY");
	mkl_set_num_threads(mkl_get_max_threads());
	mkl_disable_fast_mm();
	
	bool load_bin = false;

	SimulationData sim_data(64, 512, 64); 
	Potential pot_data(sim_data);
	WaveFunction psi(sim_data, pot_data.harmonic_trap);

	psi.get_abs(sim_data.get_N());
	psi.get_norm(sim_data);
	
	pot_data.calculate_non_linear_energy(sim_data, psi);
	pot_data.assign_position_time_evolution(sim_data, psi, true, false);
	pot_data.assign_momentum_time_evolution(sim_data, psi, false);

	if (load_bin) {
		load_binary(sim_data, psi, "ground_state.bin");	
	} 
	else {
		solve_imag(sim_data, psi, pot_data);
		save_binary(sim_data, psi, "ground_state.bin");
	}

	
	psi.get_abs(sim_data.get_N());
	psi.get_norm(sim_data);

	solve_ramp_trap(sim_data, psi, pot_data);

	psi.create_superposition(sim_data);

	psi.get_abs(sim_data.get_N());
	psi.get_norm(sim_data);

	pot_data.assign_position_time_evolution(sim_data, psi, false, true);
	pot_data.assign_momentum_time_evolution(sim_data, psi, true);

	solve_real(sim_data, psi, pot_data);
	
	printf("Solving complete\n");
	return 0;
}
