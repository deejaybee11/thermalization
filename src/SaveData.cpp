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

#include "mkl.h"

#include "../include/SaveData.hpp"

void save_2d_image(SimulationData &sim_data, WaveFunction &psi, const char * fits_file_name) {

	double *save_data;
	save_data = (double*)mkl_malloc(sim_data.get_num_x() * sim_data.get_num_y() * sizeof(double), 64);
	fitsfile *fptr;
	int status = 0;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {sim_data.get_num_y(), sim_data.get_num_x()};

	for (int i = 0; i < sim_data.get_num_x(); ++i) {
		for (int j = 0; j < sim_data.get_num_y(); ++j) {
			save_data[i * sim_data.get_num_y() + j] = 0;
			for (int k = 0; k < sim_data.get_num_z(); ++k) {
				save_data[i * sim_data.get_num_y() + j] += psi.abs_psi[i * sim_data.get_num_y() * sim_data.get_num_z() + j * sim_data.get_num_z() + k];
			}
		}
	}

	fits_create_file(&fptr, fits_file_name, &status);
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
	nelements = naxes[0] * naxes[1];
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, save_data, &status);
	fits_close_file(fptr, &status);
	fits_report_error(stderr, status);	


}
