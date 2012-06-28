/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 * Compile: gcc -O2 -o le2d *.c -lm
 */

#include <stdio.h>

#include "le_core.h"

#ifdef SAVE_EVERY_STEPS
const int save_every_step = 1;
#else
const int save_every_step = 0;
#endif

int main(int argc, char *argv[])
{
	if (argc != 4) {
		printf("Usage: %s nx ny steps.\n", argv[0]);
		return 1;
	}
	int i;
	le_point2 n = {atoi(argv[1]), atoi(argv[2])};
	int steps = atoi(argv[3]);
	le_task task;
	le_material mat;
	le_vec2 h = {1.0, 1.0};
	real dt = 0.3;
	le_vec2 center = {n.x / 2, n.y / 2};
	char name[1000];
	
	
	double t;
	unsigned long cc;
	
	/*
	 * Init material.
	 */
	le_init_material(2.0, 1.0, 1.5, &mat);
	
	/*
	 * Init task.
	 */
	le_init_task(&task, dt, h, mat, n, ST_AOS);
	
	/*
	 * Initial condition.
	 */
	le_set_ball(&task, center, 10.0, 1.0);
	
	cc = getCC();
	t = timer();
	for (i = 0; i < steps; i++) {
		if (save_every_step) {
			sprintf(name, "out-%06d.vtk", i);
			le_save_task(&task, name);
		}
		le_step(&task);
	}
	t = timer() - t;
	cc = getCC() - cc;
	/*printf("Total time, s:\t\t%f\n", t);
	printf("Time per step, s:\t%f\n", t / steps);*/
	printf("%d %d %d %f %ld\n", n.x, n.y, steps, t, cc);
	
	/*
	 * Save last step.
	 */
	le_save_task(&task, "result.vtk");
	
	/*
	 * Free memory.
	 */
	le_free_task(&task);
	return 0;
}
