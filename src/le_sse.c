#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#include "le_core.h"

#define TVD2_EPS 1e-6


#define ind(i, j) ((i) + (j) * t->n.x)
#define gind(i, j) (((le_node*)t->grid)[ind((i), (j))])

#define soa_ind(i, j, k) (t->grid[k * t->n.x * t->n.y + ind((i), (j))])
#define soa_vx(i, j) soa_ind(i, j, 0)
#define soa_vy(i, j) soa_ind(i, j, 1)
#define soa_sxx(i, j) soa_ind(i, j, 2)
#define soa_sxy(i, j) soa_ind(i, j, 3)
#define soa_syy(i, j) soa_ind(i, j, 4)

#define y_ind(i, j) ((i) + (j) * nx)
#define y_sse_ind(i, j, k) (grid[k * nx * t->n.y + y_ind((i), (j))])
#define y_sse_vx(i, j) y_sse_ind(i, j, 0)
#define y_sse_vy(i, j) y_sse_ind(i, j, 1)
#define y_sse_sxx(i, j) y_sse_ind(i, j, 2)
#define y_sse_sxy(i, j) y_sse_ind(i, j, 3)
#define y_sse_syy(i, j) y_sse_ind(i, j, 4)

inline sse_t s_min(sse_t a, sse_t b) { return sse_min(a, b); }
inline sse_t s_max(sse_t a, sse_t b) { return sse_max(a, b); }
inline sse_t s_max3(sse_t a, sse_t b, sse_t c) { return s_max(a, s_max(b, c)); }


#define limiter_minmod(r) (s_max(sse_set(0.0), s_min(sse_set(1.0), (r))))
#define limiter_cir(r) (sse_set(0.0))
#define limiter_superbee(r) (s_max3(sse_set(0.0), s_min(sse_set(1.0), sse_mul(sse_set(2.0), r)), s_min(sse_set(2.0), r)))

#define limiter limiter_superbee


inline sse_t tvd2_sse(const sse_t c, const sse_t u_2, const sse_t u_1, const sse_t u, const sse_t u1)
{
	const sse_t eps = sse_set(TVD2_EPS);
	sse_t r1 = sse_sub(u, u_1);
	sse_t r2 = sse_sub(u1, u);
	r1 = sse_add(r1, eps);
	r2 = sse_add(r2, eps);
	const sse_t r = sse_div(r1, r2);
	
	r1 = sse_sub(u_1, u_2);
	sse_t r2_1 = sse_sub(u, u_1);
	r1 = sse_add(r1, eps);
	r2_1 = sse_add(r2_1, eps);
	const sse_t r_1 = sse_div(r1, r2_1);
	
	const sse_t k = sse_mul(sse_set(0.5), sse_sub(sse_set(1.0), c));

	r2 = sse_mul(r2, k);
	r2_1 = sse_mul(r2_1, k);
	sse_t f12 = sse_add(u, sse_mul(limiter(r), r2));
	sse_t f_12 = sse_add(u_1, sse_mul(limiter(r_1), r2_1));
	return sse_mul(c, sse_sub(f_12, f12));
}

void le_sse_step_x(le_task *t)
{
	assert(t->stype == ST_SOA);
	assert(t->n.x % REAL_PER_SSE == 0);
	
	int i, j;
	int nx = t->n.x / REAL_PER_SSE;
	int ny = t->n.y;
	const sse_t k1 = sse_set(t->dt * t->mat.c1 / t->h.x);
	const sse_t k2 = sse_set(t->dt * t->mat.c2 / t->h.x);
	const sse_t irhoc1 = sse_set(t->mat.irhoc1);
	const sse_t irhoc2 = sse_set(t->mat.irhoc2);
	const sse_t rhoc1 = sse_set(t->mat.rhoc1);
	const sse_t rhoc2 = sse_set(t->mat.rhoc2);
	const sse_t rhoc3 = sse_set(t->mat.rhoc3);
	const sse_t half = sse_set(0.5f);
	sse_t *grid = (sse_t*)t->grid;
	
	real *w1[5], *w2[5], *w3[5], *w4[5];
	sse_t *sw1[5], *sw2[5], *sw3[5], *sw4[5];
#define x_malloc(w)\
	w[0] = (real*)malloc(sizeof(real) * t->n.x);\
	w[1] = (real*)malloc(sizeof(real) * t->n.x);\
	w[2] = (real*)malloc(sizeof(real) * t->n.x);\
	w[3] = (real*)malloc(sizeof(real) * t->n.x);\
	w[4] = (real*)malloc(sizeof(real) * t->n.x);\
	s##w[0] = (sse_t*)w[0];\
	s##w[1] = (sse_t*)w[1];\
	s##w[2] = (sse_t*)w[2];\
	s##w[3] = (sse_t*)w[3];\
	s##w[4] = (sse_t*)w[4];\
	
	x_malloc(w1);
	x_malloc(w2);
	x_malloc(w3);
	x_malloc(w4);
#undef x_malloc

	
#define soa_omega_x(i, j, k) \
	{ \
	const sse_t nv = y_sse_vx(i, j); \
	const sse_t N00T = sse_mul(y_sse_sxx(i, j), irhoc1); \
	const sse_t n1v = y_sse_vy(i, j); \
	const sse_t N01T = sse_mul(y_sse_sxy(i, j), irhoc2); \
	\
	sw1[k + 2][i] = sse_sub(nv, N00T); \
	sw2[k + 2][i] = sse_add(nv, N00T); \
	sw3[k + 2][i] = sse_sub(n1v, N01T); \
	sw4[k + 2][i] = sse_add(n1v, N01T); \
	}
	
	

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			soa_omega_x(i, j, 0);
		}
#define x_init(w)\
		for (i = 1; i < t->n.x; i++) {\
			w[3][i - 1] = w[2][i];\
		}\
		w[3][t->n.x - 1] = w[2][t->n.x - 1];\
		\
		for (i = 0; i < t->n.x - 1; i++) {\
			w[1][i + 1] = w[2][i];\
		}\
		w[1][0] = w[2][0];\
		for (i = 0; i < t->n.x - 2; i++) {\
			w[0][i + 2] = w[2][i];\
		}\
		w[0][0] = w[0][1] = w[2][0];\
		for (i = 2; i < t->n.x; i++) {\
			w[4][i - 2] = w[2][i];\
		}\
		w[4][t->n.x - 1] = w[4][t->n.x - 2] = w[2][t->n.x - 1];\
	
		x_init(w1);
		x_init(w2);
		x_init(w3);
		x_init(w4);
#undef x_init

		for (i = 0; i < nx; i++) {
			sse_t d1 = tvd2_sse(k1, sw1[0][i], sw1[1][i], sw1[2][i], sw1[3][i]);
			sse_t d2 = tvd2_sse(k1, sw2[4][i], sw2[3][i], sw2[2][i], sw2[1][i]);
			sse_t d3 = tvd2_sse(k2, sw3[0][i], sw3[1][i], sw3[2][i], sw3[3][i]);
			sse_t d4 = tvd2_sse(k2, sw4[4][i], sw4[3][i], sw4[2][i], sw4[1][i]);
			d1 = sse_mul(d1, half);
			d2 = sse_mul(d2, half);
			d3 = sse_mul(d3, half);
			d4 = sse_mul(d4, half);

			y_sse_vx(i, j) = sse_add(y_sse_vx(i, j), sse_add(d1, d2));
			y_sse_vy(i, j) = sse_add(y_sse_vy(i, j), sse_add(d3, d4));

			y_sse_sxx(i, j) = sse_add(y_sse_sxx(i, j), sse_mul(sse_sub(d2, d1), rhoc1));
			y_sse_syy(i, j) = sse_add(y_sse_syy(i, j), sse_mul(sse_sub(d2, d1), rhoc3));
			y_sse_sxy(i, j) = sse_add(y_sse_sxy(i, j), sse_mul(sse_sub(d4, d3), rhoc2));
		}
	}
#define x_free(w)\
	free(w[0]);\
	free(w[1]);\
	free(w[2]);\
	free(w[3]);\
	free(w[4]);
	
	x_free(w1);
	x_free(w2);
	x_free(w3);
	x_free(w4);
#undef x_free

}


void le_sse_step_y(le_task *t)
{
	assert(t->stype == ST_SOA);
	assert(t->n.x % REAL_PER_SSE == 0);
	
	int i, j;
	int nx = t->n.x / REAL_PER_SSE;
	int ny = t->n.y;
	const sse_t k1 = sse_set(t->dt * t->mat.c1 / t->h.y);
	const sse_t k2 = sse_set(t->dt * t->mat.c2 / t->h.y);
	const sse_t irhoc1 = sse_set(t->mat.irhoc1);
	const sse_t irhoc2 = sse_set(t->mat.irhoc2);
	const sse_t rhoc1 = sse_set(t->mat.rhoc1);
	const sse_t rhoc2 = sse_set(t->mat.rhoc2);
	const sse_t rhoc3 = sse_set(t->mat.rhoc3);
	const sse_t half = sse_set(0.5f);
	sse_t *grid = (sse_t*)t->grid;
	
	sse_t *w1[5], *w2[5], *w3[5], *w4[5];
#define w_malloc(w)\
	w[0] = (sse_t*)malloc(sizeof(real) * t->n.x);\
	w[1] = (sse_t*)malloc(sizeof(real) * t->n.x);\
	w[2] = (sse_t*)malloc(sizeof(real) * t->n.x);\
	w[3] = (sse_t*)malloc(sizeof(real) * t->n.x);\
	w[4] = (sse_t*)malloc(sizeof(real) * t->n.x);
	
	w_malloc(w1);
	w_malloc(w2);
	w_malloc(w3);
	w_malloc(w4);
#undef w_malloc

#define soa_omega_y(i, j, k) \
	{ \
	const sse_t nv = y_sse_vy(i, j); \
	const sse_t N00T = sse_mul(y_sse_syy(i, j), irhoc1); \
	const sse_t n1v = y_sse_vx(i, j); \
	const sse_t N01T = sse_mul(y_sse_sxy(i, j), irhoc2); \
	\
	w1[k + 2][i] = sse_sub(nv, N00T); \
	w2[k + 2][i] = sse_add(nv, N00T); \
	w3[k + 2][i] = sse_sub(n1v, N01T); \
	w4[k + 2][i] = sse_add(n1v, N01T); \
	}

	for (i = 0; i < nx; i++) {
		soa_omega_y(i, 0, 0);
		soa_omega_y(i, 1, 1);
		soa_omega_y(i, 2, 2);
	}
	
#define w_init(w)\
	for (i = 0; i < nx ; i++) {\
		w[0][i] = w[1][i] = w[2][i];\
	}
	
	w_init(w1);
	w_init(w2);
	w_init(w3);
	w_init(w4);
#undef w_init
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			sse_t d1 = tvd2_sse(k1, w1[0][i], w1[1][i], w1[2][i], w1[3][i]);
			sse_t d2 = tvd2_sse(k1, w2[4][i], w2[3][i], w2[2][i], w2[1][i]);
			sse_t d3 = tvd2_sse(k2, w3[0][i], w3[1][i], w3[2][i], w3[3][i]);
			sse_t d4 = tvd2_sse(k2, w4[4][i], w4[3][i], w4[2][i], w4[1][i]);
			d1 = sse_mul(d1, half);
			d2 = sse_mul(d2, half);
			d3 = sse_mul(d3, half);
			d4 = sse_mul(d4, half);

			y_sse_vy(i, j) = sse_add(y_sse_vy(i, j), sse_add(d1, d2));
			y_sse_vx(i, j) = sse_add(y_sse_vx(i, j), sse_add(d3, d4));

			y_sse_syy(i, j) = sse_add(y_sse_syy(i, j), sse_mul(sse_sub(d2, d1), rhoc1));
			y_sse_sxx(i, j) = sse_add(y_sse_sxx(i, j), sse_mul(sse_sub(d2, d1), rhoc3));
			y_sse_sxy(i, j) = sse_add(y_sse_sxy(i, j), sse_mul(sse_sub(d4, d3), rhoc2));
		}
		
#define w_copy(w)\
		{\
		sse_t *t = w[0];\
		w[0] = w[1];\
		w[1] = w[2];\
		w[2] = w[3];\
		w[3] = w[4];\
		w[4] = t;\
		}
		
		w_copy(w1);
		w_copy(w2);
		w_copy(w3);
		w_copy(w4);
#undef w_copy
		
		if (j < ny - 3) {
			for (i = 0; i < nx; i++) {
				soa_omega_y(i, j + 3, 2);
			}
		}
	}

#define w_free(w)\
	free(w[0]);\
	free(w[1]);\
	free(w[2]);\
	free(w[3]);\
	free(w[4]);
	
	w_free(w1);
	w_free(w2);
	w_free(w3);
	w_free(w4);
#undef w_free
}

