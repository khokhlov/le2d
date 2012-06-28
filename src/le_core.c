#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "le_core.h"


/*
 * Mapping 2d array to 1d array.
 */
#define ind(i, j) ((i) + (j) * t->n.x)
#define gind(i, j) (((le_node*)t->grid)[ind((i), (j))])

/*
 * Vector norm.
 */
#define vnorm(v) (sqrt(v.x * v.x + v.y * v.y))

#define TVD2_EPS 1e-6

inline real le_min(real a, real b) { return a > b ? b : a; }
inline real le_max(real a, real b) { return a > b ? a : b; }
inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }


#define limiter_minmod(r) (le_max(0.0, le_min(1.0, (r))))
#define limiter_cir(r) (0.0)
#define limiter_superbee(r) (le_max3(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r)))

/*
 * Set TVD limiter (http://en.wikipedia.org/wiki/Flux_limiters).
 */
#define limiter limiter_superbee

/*
 * Second order TVD scheme (http://en.wikipedia.org/wiki/Total_variation_diminishing).
 */
inline real tvd2(const real c, const real u_2, const real u_1, const real u, const real u1)
{
	real r1 = (u  - u_1);
	real r2 = (u1 - u);
	if (r2 == 0.0) {
		r1 += TVD2_EPS;
		r2 += TVD2_EPS;
	}
	const real r = r1 / r2;
	r1 = (u_1 - u_2);
	r2 = (u   - u_1);
	if (r2 == 0.0) {
		r1 += TVD2_EPS;
		r2 += TVD2_EPS;
	}
	const real r_1 = r1 / r2;

	const real f12  = u   + limiter(r)   / 2.0 * (1.0 - c) * (u1 - u);
	const real f_12 = u_1 + limiter(r_1) / 2.0 * (1.0 - c) * (u  - u_1);

	return u - c * (f12 - f_12);
}

void le_set_ball(le_task *t, const le_vec2 c, const real r, const real s)
{
	assert(t->stype == ST_AOS);
	int i, j;
	for (i = 0; i < t->n.x; i++) {
		for (j = 0; j < t->n.y; j++) {
			le_vec2 x = {t->h.x * i, t->h.y * j};
			le_vec2 d = {x.x - c.x, x.y - c.y};
			if (vnorm(d) < r) {
				/*
				 * Set pressure disturbance,
				 */
				gind(i, j).s.xx = s;
				gind(i, j).s.yy = s;
			}
		}
	}
}

/*
 * Write float to file and reverse byte order.
 */
void write_float(FILE* f, const float v)
{
	union {
		float f;
		unsigned char b[4];
	} dat1, dat2;
	dat1.f = v;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	fwrite(dat2.b, sizeof(unsigned char), 4, f);
}


void le_init_task(le_task *task, const real dt, const le_vec2 h, const le_material mat, const le_point2 n, const int stype)
{
	task->dt  = dt;
	task->h   = h;
	task->mat = mat;
	task->n   = n;
	task->grid = (real*)malloc(sizeof(le_node) * n.x * n.y);
	task->stype = stype;
	memset(task->grid, 0, sizeof(le_node) * n.x * n.y);
}

void le_free_task(le_task* task)
{
	free(task->grid);
}

int le_save_task(le_task *t, const char *file)
{
	assert(t->stype == ST_AOS);
	int i, j;
	FILE *fp = fopen(file, "w");
	if (fp == NULL) {
		perror("Failed to open file");
		return 1;
	}
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Created by le_save_task\n");
	fprintf(fp, "BINARY\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS\n");
	fprintf(fp, "DIMENSIONS %d %d 1\n", t->n.x, t->n.y);
	fprintf(fp, "SPACING %f %f 0.0\n", t->h.x, t->h.y);
	fprintf(fp, "ORIGIN 0.0 0.0 0.0\n");
	fprintf(fp, "POINT_DATA %d\n", t->n.x * t->n.y);
	
	/* velocity */
	fprintf(fp, "SCALARS v float 1\n");
	fprintf(fp, "LOOKUP_TABLE v_table\n");
	for (j = 0; j < t->n.y; j++) {
		for (i = 0; i < t->n.x; i++) {
			float v = vnorm(gind(i, j).v);
			write_float(fp, v);
		}
	}
	/*
	 * You can use the same code for saving other variables.
	 */
	fclose(fp);
	return 0;
}

void le_init_material(const real c1, const real c2, const real rho, le_material *m)
{
	m->c1 = c1;
	m->c2 = c2;
	m->rho = rho;

	/*
	 * Cached values.
	 */
	m->irhoc1 = 1.0 / (c1 * rho);
	m->irhoc2 = 1.0 / (c2 * rho);
	m->rhoc1 = c1 * rho;
	m->rhoc2 = c2 * rho;
	real mu = rho * c2 * c2;
	real la = rho * c1 * c1 - 2.0 * mu;
	m->rhoc3 = rho * c1 * la / (la + 2.0 * mu);
}

/*
 * w = OmegaR * u.
 */
inline void omega_x(const le_material *m, const le_node *n, le_w *w)
{
	const real nv = n->v.x;
	const real N00T = n->s.xx * m->irhoc1;

	const real n1v = n->v.y;
	const real N01T = n->s.xy * m->irhoc2;

	w->w1 = nv  - N00T;
	w->w2 = nv  + N00T;
	w->w3 = n1v - N01T;
	w->w4 = n1v + N01T;
}

inline void omega_y(const le_material *m, const le_node *n, le_w *w)
{
	const real nv = n->v.y;
	const real N00T = n->s.yy * m->irhoc1;

	const real n1v = n->v.x;
	const real N01T = n->s.xy * m->irhoc2;

	w->w1 = nv  - N00T;
	w->w2 = nv  + N00T;
	w->w3 = n1v - N01T;
	w->w4 = n1v + N01T;
}


inline void inc_x(const le_material *m, le_node *n, const le_w *d)
{
	const real d1 = 0.5 * d->w1;
	const real d2 = 0.5 * d->w2;
	const real d3 = 0.5 * d->w3;
	const real d4 = 0.5 * d->w4;

	n->v.x += d1 + d2;
	n->v.y += d3 + d4;

	n->s.xx += (d2 - d1) * m->rhoc1;
	n->s.yy += (d2 - d1) * m->rhoc3;
	n->s.xy += m->rhoc2 * (d4 - d3);
}

inline void inc_y(const le_material *m, le_node *n, const le_w *d)
{
	const real d1 = 0.5 * d->w1;
	const real d2 = 0.5 * d->w2;
	const real d3 = 0.5 * d->w3;
	const real d4 = 0.5 * d->w4;

	n->v.y += d1 + d2;
	n->v.x += d3 + d4;

	n->s.yy += (d2 - d1) * m->rhoc1;
	n->s.xx += (d2 - d1) * m->rhoc3;
	n->s.xy += m->rhoc2 * (d4 - d3);
}

inline void reconstruct(const le_w ppu, const le_w pu, const le_w u, const le_w nu, const le_w nnu, const real k1, const real k2, le_w *d)
{
	d->w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1) - u.w1; // c1
	d->w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2) - u.w2; // -c1
	d->w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3) - u.w3; // c2
	d->w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4) - u.w4; // -c2
}

void le_step_x(le_task *t)
{
	assert(t->stype == ST_AOS);
	/*
	 * Due to our system of pde is linear, we can use some simple way to solve it.
	 * du/dt + A * du/dx = 0.
	 * Vector u = {vx, vy, sxx, sxy, syy}.
	 * Matrix A could be represent in form OmegaL * Lambda * OmegaR,
	 * where Lambda - diagonal matrix of eigen values of matrix A.
	 * In our case Lambda = diag{c1, -c1, c2, -c2, 0}.
	 * OmegaR and OmegaL - metrices from eigen vectors of matrix A,
	 * OmegaR * OmegaL = E, where E = diag{1, 1, 1, 1, 1}.
	 * 
	 * We can rewrite out system in form:
	 * du/dt + OmegaL * Lambda * OmegaR du/dx = 0, multiply on matrix OmegaR:
	 * 
	 * OmegaR * du/dt + OmegaR * OmegaL * Lambda * OmegaR du/dx = 0.
	 * 
	 * Introduce new variables (http://en.wikipedia.org/wiki/Riemann_invariant):
	 * w = {w1, w2, w3, w4, w5},
	 * w = OmegaR * u, then we got:
	 * 
	 * dw/dt + Lambda * dw/dx = 0.
	 * And we get system of independent advection equations, that we can solve separatly.
	 * 
	 * So we get next algorithm:
	 * 1. Introduce new variables w = OmegaR * u;
	 * 2. Solve 5 equations of linear advection (in real we solve only 4, because of in fifth equation speed is 0);
	 * 3. Make inverse transformation u = OmegaL * w. 
	 */
	int i, j;
	
	// Courant number (http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition).
	const real k1 = t->dt * t->mat.c1 / t->h.x;
	const real k2 = t->dt * t->mat.c2 / t->h.x;
	
	for (j = 0; j < t->n.y; j++) {
		/*
		 * Riemann invariants for 5-point sctencil difference scheme.
		 */
		le_w w_2, w_1, w, w1, w2;
		
		omega_x(&t->mat, &gind(0, j), &w);
		omega_x(&t->mat, &gind(1, j), &w1);
		omega_x(&t->mat, &gind(2, j), &w2);
		
		w_2 = w_1 = w; // In linear case we can do this.
		
		for (i = 0; i < t->n.x; i++) {
			le_w d;
			reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
			inc_x(&t->mat, &gind(i, j), &d);
			w_2 = w_1;
			w_1 = w;
			w   = w1;
			w1  = w2;
			if (i < t->n.x - 3) omega_x(&t->mat, &gind(i + 3, j), &w2);
		}
	}
	
}

void le_step_y(le_task *t)
{
	assert(t->stype == ST_AOS);
	/*
	 * The same way as step_x.
	 * WARNING: Not cache friendly!
	 */
	int i, j;
	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;
	
	for (i = 0; i < t->n.x; i++) {
		le_w w_2, w_1, w, w1, w2;
		
		omega_y(&t->mat, &gind(i, 0), &w);
		omega_y(&t->mat, &gind(i, 1), &w1);
		omega_y(&t->mat, &gind(i, 2), &w2);
		
		w_2 = w_1 = w;
		
		for (j = 0; j < t->n.y; j++) {
			le_w d;
			reconstruct(w_2, w_1, w, w1, w2, k1, k2, &d);
			inc_y(&t->mat, &gind(i, j), &d);
			w_2 = w_1;
			w_1 = w;
			w   = w1;
			w1  = w2;
			if (j < t->n.y - 3) omega_y(&t->mat, &gind(i, j + 3), &w2);
		}
	}
}

void le_step_y_cf(le_task *t, const int_t cfs)
{
	assert(t->stype == ST_AOS);
	assert(cfs > 0 && cfs <= t->n.x);
	
	/*
	 * Cache friendly version of function le_step_y.
	 */
	int i, j, k;
	const real k1 = t->dt * t->mat.c1 / t->h.y;
	const real k2 = t->dt * t->mat.c2 / t->h.y;
	
	le_w *w_2, *w_1, *w, *w1, *w2;
	w_2 = (le_w*)malloc(sizeof(le_w) * cfs);
	w_1 = (le_w*)malloc(sizeof(le_w) * cfs);
	w   = (le_w*)malloc(sizeof(le_w) * cfs);
	w1  = (le_w*)malloc(sizeof(le_w) * cfs);
	w2  = (le_w*)malloc(sizeof(le_w) * cfs);
	

	for (k = 0; k < t->n.x; k += cfs) {
		int cfs_n = (k + cfs > t->n.x) ? t->n.x - k : cfs;
		/* Prepare values and copy it to w[] */
		for (i = 0; i < cfs_n; i++) {
			omega_y(&t->mat, &gind(i + k, 0), w  + i);
			omega_y(&t->mat, &gind(i + k, 1), w1 + i);
			omega_y(&t->mat, &gind(i + k, 2), w2 + i);
		}
		memcpy(w_2, w, sizeof(le_w) * cfs_n);
		memcpy(w_1, w, sizeof(le_w) * cfs_n);
		
		/*
		 * Iterate throw cached arrays of values.
		 */
		for (j = 0; j < t->n.y; j++) {
			for (i = 0; i < cfs_n; i++) {
				le_w d;
				reconstruct(w_2[i], w_1[i], w[i], w1[i], w2[i], k1, k2, &d);
				inc_y(&t->mat, &gind(i + k, j), &d);
			}
			/*
			 * Swap values.
			 */
			le_w *wt = w_2;
			w_2 = w_1;
			w_1 = w;
			w   = w1;
			w1  = w2;
			w2  = wt;
			if (j < t->n.y - 3) {
				for (i = 0; i < cfs_n; i++) {
					omega_y(&t->mat, &gind(i + k, j + 3), w2 + i);
				}
			}
		}
	}
	
	free(w_2);
	free(w_1);
	free(w);
	free(w1);
	free(w2);
}


void le_step(le_task *task)
{
	/*
	 * We solve regular hyperbolic system of PDE (http://en.wikipedia.org/wiki/Hyperbolic_partial_differential_equation) in form:
	 * du/dt + Ax * du/dx + Ay * du/dy = 0.
	 * 
	 * During time integration we use dimension split method:
	 * 1. Step:
	 * Integrate system dv/dt + Ax * dv/dx = 0, get u = v^(n + 1).
	 * 2. Step:
	 * Integrate system du/dt + Ay * du/dy = 0, get on next time step u^(n + 1).
	 */
	
	le_step_x(task);
	le_step_y(task);
}

void le_step_cf(le_task *task, const int cfs)
{
	le_step_x(task);
	le_step_y_cf(task, cfs);
}


