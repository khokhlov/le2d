/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012
 */

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

#include "le_core.h"

#ifdef SAVE_EVERY_STEPS
const int save_every_step = 1;
#else
const int save_every_step = 0;
#endif

void le_sse_step_y(le_task *t);
void le_sse_step_x(le_task *t);
void le_step_sse(le_task *t)
{
	le_sse_step_x(t);
	le_sse_step_y(t);
}

typedef struct {
	int rank;
	int steps;
	pthread_barrier_t *b;
	pthread_t th;
	le_task task;
} le_pthread_task;

void *tfunc(void *arg)
{
	le_pthread_task *t = (le_pthread_task*)arg;
	int i;
	char name[100];
	for (i = 0; i < t->steps; i++) {
		pthread_barrier_wait(t->b);
		if (save_every_step) {
			if (t->rank == 0) {
				sprintf(name, "out-%06d.vtk", i);
				le_save_task(&t->task, name);
			}
		}
		le_step_sse(&t->task);
	}
}

int main(int argc, char *argv[])
{
	if (argc != 5) {
		printf("Usage: %s nx ny steps num_threads.\n", argv[0]);
		return 1;
	}
	int i;
	le_point2 n = {atoi(argv[1]), atoi(argv[2])};
	int steps = atoi(argv[3]);
	int tnum = atoi(argv[4]);
	le_task task;
	le_material mat;
	le_vec2 h = {1.0, 1.0};
	real dt = 0.3;
	le_vec2 center = {n.x / 2, n.y / 2};
	char name[1000];
	pthread_barrier_t bar;
	
	/*
	 * Create threads array.
	 */
	le_pthread_task *th = (le_pthread_task*)malloc(sizeof(le_pthread_task) * tnum);
	pthread_barrier_init(&bar, NULL, tnum);
	
	
	double t;
	unsigned long cc;
	
	/*
	 * Init material.
	 */
	le_init_material(2.0, 1.0, 1.5, &mat);
	
	/*
	 * Init task.
	 */
	le_init_task(&task, dt, h, mat, n, ST_SOA);
	
	/*
	 * Initial condition.
	 */
	le_set_ball(&task, center, 10.0, 1.0);
	
	
	
	cc = getCC();
	t = timer();
	for (i = 0; i < tnum; i++) {
		th[i].task = task;
		th[i].task.jmin = (task.n.y / tnum) * i + (i > task.n.y % tnum);
		th[i].task.jmax = (task.n.y / tnum) * (i + 1) + ((i + 1) > task.n.y % tnum);
		th[i].rank = i;
		th[i].steps = steps;
		th[i].b = &bar;
		pthread_create(&th[i].th, NULL, tfunc, &th[i]);
	}
	for (i = 0; i < tnum; i++) {
		pthread_join(th[i].th, NULL);
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
	free(th);
	pthread_barrier_destroy(&bar);
	return 0;
}
