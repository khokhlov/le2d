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

inline sse_t sse_min(sse_t a, sse_t b) { return _mm_min_pd(a, b); }
inline sse_t sse_max(sse_t a, sse_t b) { return _mm_max_pd(a, b); }
inline sse_t sse_max3(sse_t a, sse_t b, sse_t c) { return sse_max(a, sse_max(b, c)); }


#define limiter_minmod(r) (sse_max(_mm_set1_pd(0.0), sse_min(_mm_set1_pd(1.0), (r))))
#define limiter_cir(r) (_mm_set1_pd(0.0))
#define limiter_superbee(r) (sse_max3(_mm_set1_pd(0.0), sse_min(_mm_set1_pd(1.0), _mm_mul_pd(_mm_set1_pd(2.0), r)), sse_min(_mm_set1_pd(2.0), r)))

#define limiter limiter_superbee


inline sse_t tvd2_sse(const sse_t c, const sse_t u_2, const sse_t u_1, const sse_t u, const sse_t u1)
{
	const sse_t eps = _mm_set1_pd(TVD2_EPS);
	sse_t r1 = _mm_sub_pd(u, u_1);
	sse_t r2 = _mm_sub_pd(u1, u);
	r1 = _mm_add_pd(r1, eps);
	r2 = _mm_add_pd(r2, eps);
	const sse_t r = _mm_div_pd(r1, r2);
	
	r1 = _mm_sub_pd(u_1, u_2);
	sse_t r2_1 = _mm_sub_pd(u, u_1);
	r1 = _mm_add_pd(r1, eps);
	r2_1 = _mm_add_pd(r2_1, eps);
	const sse_t r_1 = _mm_div_pd(r1, r2_1);
	
	const sse_t k = _mm_mul_pd(_mm_set1_pd(0.5), _mm_sub_pd(_mm_set1_pd(1.0), c));

	r2 = _mm_mul_pd(r2, k);
	r2_1 = _mm_mul_pd(r2_1, k);
	sse_t f12 = _mm_add_pd(u, _mm_mul_pd(limiter(r), r2));
	sse_t f_12 = _mm_add_pd(u_1, _mm_mul_pd(limiter(r_1), r2_1));
	return _mm_mul_pd(c, _mm_sub_pd(f_12, f12));
}

void le_sse_step_y(le_task *t)
{
	assert(t->stype == ST_SOA);
	assert(t->n.x % REAL_PER_SSE == 0);
	
	int i, j;
	int nx = t->n.x / REAL_PER_SSE;
	int ny = t->n.y;
	const sse_t k1 = _mm_set1_pd(t->dt * t->mat.c1 / t->h.y);
	const sse_t k2 = _mm_set1_pd(t->dt * t->mat.c2 / t->h.y);
	const sse_t irhoc1 = _mm_set1_pd(t->mat.irhoc1);
	const sse_t irhoc2 = _mm_set1_pd(t->mat.irhoc2);
	const sse_t rhoc1 = _mm_set1_pd(t->mat.rhoc1);
	const sse_t rhoc2 = _mm_set1_pd(t->mat.rhoc2);
	const sse_t rhoc3 = _mm_set1_pd(t->mat.rhoc3);
	const sse_t half = _mm_set1_pd(0.5f);
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
	const sse_t N00T = _mm_mul_pd(y_sse_syy(i, j), irhoc1); \
	const sse_t n1v = y_sse_vx(i, j); \
	const sse_t N01T = _mm_mul_pd(y_sse_sxy(i, j), irhoc2); \
	\
	w1[k + 2][i] = _mm_sub_pd(nv, N00T); \
	w2[k + 2][i] = _mm_add_pd(nv, N00T); \
	w3[k + 2][i] = _mm_sub_pd(n1v, N01T); \
	w4[k + 2][i] = _mm_add_pd(n1v, N01T); \
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
			d1 = _mm_mul_pd(d1, half);
			d2 = _mm_mul_pd(d2, half);
			d3 = _mm_mul_pd(d3, half);
			d4 = _mm_mul_pd(d4, half);

			y_sse_vy(i, j) = _mm_add_pd(y_sse_vy(i, j), _mm_add_pd(d1, d2));
			y_sse_vx(i, j) = _mm_add_pd(y_sse_vx(i, j), _mm_add_pd(d3, d4));

			y_sse_syy(i, j) = _mm_add_pd(y_sse_syy(i, j), _mm_mul_pd(_mm_sub_pd(d2, d1), rhoc1));
			y_sse_sxx(i, j) = _mm_add_pd(y_sse_sxx(i, j), _mm_mul_pd(_mm_sub_pd(d2, d1), rhoc3));
			y_sse_sxy(i, j) = _mm_add_pd(y_sse_sxy(i, j), _mm_mul_pd(_mm_sub_pd(d3, d4), rhoc2));
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

