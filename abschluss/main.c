#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>

/*
Compile with: gcc -O2 -fopenmp -march=native main.c -o main
*/

/* Jacobi */
void jacobiSequential(const float* startVector, float h, const float* functionTable, float* jacobiResult);
void jacobi(const float* startVector, float h, const float* functionTable, float* jacobiResult);
void jacobiSSE(const float* startVector, float h, const float* functionTable, float* jacobiResult);

/* Gauss Seidel */
void gaussSeidel(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
/* Naive parallelization */
void gaussSeidelNaiv(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
/* RotSchwarz parallelization */
void gaussSeidelRotSchwarz(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzEven(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzOdd(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzSSE(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzEvenSSE(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzOddSSE(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
/* Wavefront parallelization */
void gaussSeidelWavefront(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelWavefrontCache(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);

/* Helper functions */
float f1(float x, float y);
void computeFunctionTable(float h, float* functionTable);
void printResultMatrix(const float* matrix);
void printAnalyticalResult(float h);
bool compare(float* m1, float* m2);
inline __m128 abs_ps(__m128 x);
inline float sse_sum(__m128 x);


int size;

#define CO(i,j) ( (j) * (size) + (i) )

const static int MAX_ITERATIONS = 999999999;
const static float EPSILON = 0.01;
const static float TOL = 0.00000001;

inline float sse_sum(__m128 x) {
    float res;
    x = _mm_hadd_ps(x, x);
    x = _mm_hadd_ps(x, x);
    _mm_store_ss(&res, x);
    return res;
}

inline __m128 abs_ps(__m128 x) {
    const __m128 sign_mask = _mm_set1_ps(-0.f); // -0.f = 1 << 31
    return _mm_andnot_ps(sign_mask, x);
}

double get_wall_time()   // returns wall time in seconds
{
    struct timeval time;
    if (gettimeofday(&time,NULL))
    {
        return 0;
    }
    double wtime = (double)time.tv_sec + (double)time.tv_usec * 0.000001;
    return wtime;
}
void jacobiSerial(const float* startVector, float h, const float* functionTable, float* jacobiResult)
{
    int i, j, k;
    float* array0 = (float*) malloc(size * size * sizeof(float));
    float* array1 = (float*) malloc(size * size * sizeof(float));

    float* a0 = array0; // last iteration
    float* a1 = array1; // current iteration


    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;

        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;

        for (j = 1; j < size - 1; j++)
        {
            for (i = 1; i < size - 1; i++)
            {
                a1[CO(i,j)] = a0[CO(i, j - 1)]
                              + a0[CO(i - 1, j)]
                              + a0[CO(i, j + 1)]
                              + a0[CO(i + 1, j)]
                              + functionTable[CO(i, j)];
                a1[CO(i,j)] *= 0.25;

                diff += fabsf(a1[CO(i,j)] - a0[CO(i,j)]);
            }
        }

        if (diff / (size * size) < TOL) break;
    }


    for (i = 0; i < size * size; ++i)
    {
        jacobiResult[i] = a1[i];
    }

    free(a0);
    free(a1);
}

void jacobi(const float* startVector, float h, const float* functionTable, float* jacobiResult)
{
    int i, j, k;
    float* array0 = (float*) malloc(size * size * sizeof(float));
    float* array1 = (float*) malloc(size * size * sizeof(float));

    float* a0 = array0; // last iteration
    float* a1 = array1; // current iteration

    #pragma omp parallel for private(i)
    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;

        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;

        #pragma omp parallel for private(j, i) collapse(2) reduction(+:diff)
        for (j = 1; j < size - 1; j++)
        {
            for (i = 1; i < size - 1; i++)
            {
                a1[CO(i,j)] = a0[CO(i, j - 1)]
                              + a0[CO(i - 1, j)]
                              + a0[CO(i, j + 1)]
                              + a0[CO(i + 1, j)]
                              + functionTable[CO(i, j)];
                a1[CO(i,j)] *= 0.25;

                diff += fabsf(a1[CO(i,j)] - a0[CO(i,j)]);
            }
        }

        if (diff / (size * size) < TOL) break;
    }

    #pragma omp parallel for private(i)
    for (i = 0; i < size * size; ++i)
    {
        jacobiResult[i] = a1[i];
    }

    free(a0);
    free(a1);
}

void jacobiSSE(const float* startVector, float h, const float* functionTable, float* jacobiResult)
{
    const __m128 vec_0_25 = _mm_set_ps1(0.25);
    int i, j, k;
    float* array0 = (float*) malloc(size * size * sizeof(float));
    float* array1 = (float*) malloc(size * size * sizeof(float));

    float* a0 = array0; // last iteration
    float* a1 = array1; // current iteration

    #pragma omp parallel for private(i)
    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;
        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;

        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j++)
        {
            for (i = 1; i < size - 5; i += 4)
            {
    		__m128 vec_left = _mm_loadu_ps((float const*) &a0[CO(i, j - 1)]);
    		__m128 vec_up = _mm_loadu_ps((float const*) &a0[CO(i - 1, j)]);
    		__m128 vec_right = _mm_loadu_ps((float const*) &a0[CO(i, j + 1)]);
    		__m128 vec_down = _mm_loadu_ps((float const*) &a0[CO(i + 1, j)]);
    		__m128 vec_ft = _mm_set_ps(functionTable[CO(i+3, j)], functionTable[CO(i+2, j)],
    				functionTable[CO(i+1, j)], functionTable[CO(i, j)]);

    		__m128 vec_a1 = _mm_add_ps(vec_left, vec_up);
    		vec_a1 = _mm_add_ps(vec_a1, vec_right);
    		vec_a1 = _mm_add_ps(vec_a1, vec_down);
    		vec_a1 = _mm_add_ps(vec_a1, vec_ft);
    		vec_a1 = _mm_mul_ps(vec_a1, vec_0_25);
    		_mm_storeu_ps(&a1[CO(i,j)], vec_a1);

    		__m128 vec_old = _mm_loadu_ps((float const*) &a0[CO(i, j)]);
    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_a1, vec_old));

    		diff += sse_sum(vec_diff);
            }

            for (; i < size - 1; i++) // do the rest sequentially
            {
                a1[CO(i,j)] = a0[CO(i, j - 1)]
                              + a0[CO(i - 1, j)]
                              + a0[CO(i, j + 1)]
                              + a0[CO(i + 1, j)]
                              + functionTable[CO(i, j)];
                a1[CO(i,j)] *= 0.25;

                diff += fabsf(a1[CO(i,j)] - a0[CO(i,j)]);
            }
        }

        if (diff / (size * size) < TOL) break;
    }

    #pragma omp parallel for private(i)
    for (i = 0; i < size * size; ++i)
    {
        jacobiResult[i] = a1[i];
    }

    free(a0);
    free(a1);
}

void gaussSeidel(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    int i, j, k;

    float* array0 = (float*) malloc(size * size * sizeof(float));
    float* array1 = (float*) malloc(size * size * sizeof(float));

    float* a0 = array0; // last iteration
    float* a1 = array1; // current iteration

    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;

        float diff = 0;
        for (j = 1; j < size - 1; j++)
        {
            for (i = 1; i < size - 1; i++)
            {
                a1[CO(i,j)] = a1[CO(i, j - 1)]
                              + a1[CO(i - 1, j)]
                              + a0[CO(i, j + 1)]
                              + a0[CO(i + 1, j)]
                              + functionTable[CO(i, j)];
                a1[CO(i,j)] *= 0.25;

		diff += fabsf(a1[CO(i,j)] - a0[CO(i,j)]);
            }
        }

        if (diff / (size * size) < TOL) break;
    }

    for (i = 0; i < size * size; ++i)
    {
        gaussSeidelResult[i] = a1[i];
    }
    free(a0);
    free(a1);
}

void gaussSeidelRotSchwarz(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    if (size % 2 == 0)
    {
        gaussSeidelRotSchwarzEven(startVector, h, functionTable, gaussSeidelResult);
    }
    else
    {
        gaussSeidelRotSchwarzOdd(startVector, h, functionTable, gaussSeidelResult);
    }
}

void gaussSeidelRotSchwarzSSE(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    if (size % 2 == 0)
    {
        gaussSeidelRotSchwarzEvenSSE(startVector, h, functionTable, gaussSeidelResult);
    }
    else
    {
        gaussSeidelRotSchwarzOddSSE(startVector, h, functionTable, gaussSeidelResult);
    }
}

// For size % 2 == 0
void gaussSeidelRotSchwarzEven(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    int i1, i, j, k;
    int halfSize = size / 2;
    float* r1 = (float*) malloc(size * halfSize * sizeof(float));
    float* s1 = (float*) malloc(size * halfSize * sizeof(float));

    // fill the arrays
    #pragma omp parallel for private(i, j)
    for (j = 0; j < size-1; j+=2)
    {
        for (i = 0; i < halfSize; ++i) // for even j
        {
            int baseIdx = j * halfSize;
            int idx = baseIdx + i;
            int idxRot, idxSchwarz;

            idxRot = 2 * idx;
            idxSchwarz = 2 * idx + 1;


            r1[idx] = startVector[idxRot];
            s1[idx] = startVector[idxSchwarz];
        }

        for (i = 0; i < halfSize; ++i) // for odd j+1
        {
            int baseIdx = (j+1) * halfSize;
            int idx = baseIdx + i;
            int idxRot, idxSchwarz;

            idxRot = 2 * idx + 1;
            idxSchwarz = 2 * idx;

            r1[idx] = startVector[idxRot];
            s1[idx] = startVector[idxSchwarz];
        }
    }

    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;
        // rote Punkte
        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 1; i < size - 2; i += 2) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		float old_val = r1[idx];

			r1[idx] = s1[idx - halfSize] // links
		                  + s1[idx] // oben
		                  + s1[idx + halfSize] // rechts
		                  + s1[idx + 1] // unten
		                  + functionTable[idxWhole];
		        r1[idx] *= 0.25;

		        diff += fabsf(r1[idx] - old_val);
	    	}

	    	for (i = 2; i < size - 1; i += 2) {
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
	    		float old_val = r1[idx];

			r1[idx] = s1[idx - halfSize] // links
			          + s1[idx - 1] // oben
			          + s1[idx + halfSize] // rechts
			          + s1[idx] // unten
			          + functionTable[idxWhole];
			r1[idx] *= 0.25;

			diff += fabsf(r1[idx] - old_val);
	    	}
        }

        // schwarze Punkte
        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 2; i < size - 1; i += 2) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		float old_val = s1[idx];

	    		s1[idx] = r1[idx - halfSize] // links
		                  + r1[idx - 1] // oben
		                  + r1[idx + halfSize] // rechts
		                  + r1[idx] // unten
		                  + functionTable[idxWhole];
		        s1[idx] *= 0.25;

		        diff += fabsf(s1[idx] - old_val);
	    	}

	    	for (i = 1; i < size - 2; i += 2) {
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
	    		float old_val = s1[idx];

			s1[idx] = r1[idx - halfSize] // links
			          + r1[idx] // oben
			          + r1[idx + halfSize] // rechts
			          + r1[idx + 1] // unten
			          + functionTable[idxWhole];
			s1[idx] *= 0.25;

			diff += fabsf(s1[idx] - old_val);
	    	}
        }

        if (diff / (size * size) < TOL) break;
    }

    #pragma omp parallel for private(i, j)
    for (j = 0; j < size-1; j+=2)
    {
        for (i = 0; i < halfSize; ++i) // for even j
        {
            int baseIdx = j * halfSize;
            int idx = baseIdx + i;
            int idxRot = 2 * idx;
            int idxSchwarz = 2 * idx + 1;
            gaussSeidelResult[idxRot] = r1[idx];
            gaussSeidelResult[idxSchwarz] = s1[idx];
        }
        for (i = 0; i < halfSize; ++i) // for odd j+1
        {
            int baseIdx = (j+1) * halfSize;
            int idx = baseIdx + i;
            int idxRot = 2 * idx + 1;
            int idxSchwarz = 2 * idx;
            gaussSeidelResult[idxRot] = r1[idx];
            gaussSeidelResult[idxSchwarz] = s1[idx];
        }
    }

    free(r1);
    free(s1);
}

// For size % 2 == 0
void gaussSeidelRotSchwarzEvenSSE(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    const __m128 vec_0_25 = _mm_set_ps1(0.25);
    int i1, i, j, k;
    int halfSize = size / 2;
    float* r1 = (float*) malloc(size * halfSize * sizeof(float));
    float* s1 = (float*) malloc(size * halfSize * sizeof(float));

    // fill the arrays
    #pragma omp parallel for private(i, j)
    for (j = 0; j < size-1; j+=2)
    {
        for (i = 0; i < halfSize; ++i) // for even j
        {
            int baseIdx = j * halfSize;
            int idx = baseIdx + i;
            int idxRot, idxSchwarz;

            idxRot = 2 * idx;
            idxSchwarz = 2 * idx + 1;


            r1[idx] = startVector[idxRot];
            s1[idx] = startVector[idxSchwarz];
        }

        for (i = 0; i < halfSize; ++i) // for odd j+1
        {
            int baseIdx = (j+1) * halfSize;
            int idx = baseIdx + i;
            int idxRot, idxSchwarz;

            idxRot = 2 * idx + 1;
            idxSchwarz = 2 * idx;

            r1[idx] = startVector[idxRot];
            s1[idx] = startVector[idxSchwarz];
        }
    }

    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;
        // rote Punkte
        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j+=2)
        {
        	for (i = 1; i < size - 8; i += 8) {
        		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		__m128 vec_old = _mm_loadu_ps((float const*) &r1[idx]);

	    		__m128 vec_left = _mm_loadu_ps((float const*) &s1[idx - halfSize]);
	    		__m128 vec_up = _mm_loadu_ps((float const*) &s1[idx]);
	    		__m128 vec_right = _mm_loadu_ps((float const*) &s1[idx + halfSize]);
	    		__m128 vec_down = _mm_loadu_ps((float const*) &s1[idx + 1]);
	    		__m128 vec_ft = _mm_set_ps(functionTable[idxWhole + 6], functionTable[idxWhole + 4],
	    				functionTable[idxWhole + 2], functionTable[idxWhole]);

	    		__m128 vec_r1 = _mm_add_ps(vec_left, vec_up);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_right);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_down);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_ft);
	    		vec_r1 = _mm_mul_ps(vec_r1, vec_0_25);
	    		_mm_storeu_ps(&r1[idx], vec_r1);

	    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_r1, vec_old));
	    		diff += sse_sum(vec_diff);
	    	}

	    	for (; i < size - 2; i += 2) { // do the remainder sequentially
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		const float old_val = r1[idx];

			r1[idx] = s1[idx - halfSize] // links
		                  + s1[idx] // oben
		                  + s1[idx + halfSize] // rechts
		                  + s1[idx + 1] // unten
		                  + functionTable[idxWhole];
		        r1[idx] *= 0.25;

		        diff += fabsf(r1[idx] - old_val);
	    	}

	    	for (i = 2; i < size - 7; i += 8) {
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
	    		__m128 vec_old = _mm_loadu_ps((float const*) &r1[idx]);

	    		__m128 vec_left = _mm_loadu_ps((float const*) &s1[idx - halfSize]);
	    		__m128 vec_up = _mm_loadu_ps((float const*) &s1[idx - 1]);
	    		__m128 vec_right = _mm_loadu_ps((float const*) &s1[idx + halfSize]);
	    		__m128 vec_down = _mm_loadu_ps((float const*) &s1[idx]);
	    		__m128 vec_ft = _mm_set_ps(functionTable[idxWhole + 6], functionTable[idxWhole + 4],
	    				functionTable[idxWhole + 2], functionTable[idxWhole]);

	    		__m128 vec_r1 = _mm_add_ps(vec_left, vec_up);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_right);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_down);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_ft);
	    		vec_r1 = _mm_mul_ps(vec_r1, vec_0_25);
	    		_mm_storeu_ps(&r1[idx], vec_r1);

	    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_r1, vec_old));
	    		diff += sse_sum(vec_diff);
	    	}

	    	for (i = 2; i < size - 1; i += 2) { // do the remainder sequentially
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
	    		float old_val = r1[idx];

			r1[idx] = s1[idx - halfSize] // links
			          + s1[idx - 1] // oben
			          + s1[idx + halfSize] // rechts
			          + s1[idx] // unten
			          + functionTable[idxWhole];
			r1[idx] *= 0.25;
			diff += fabsf(r1[idx] - old_val);
	    	}
        }

        // schwarze Punkte
        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 2; i < size - 7; i += 8) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		__m128 vec_old = _mm_loadu_ps((float const*) &s1[idx]);

		        __m128 vec_left = _mm_loadu_ps((float const*) &r1[idx - halfSize]);
	    		__m128 vec_up = _mm_loadu_ps((float const*) &r1[idx - 1]);
	    		__m128 vec_right = _mm_loadu_ps((float const*) &r1[idx + halfSize]);
	    		__m128 vec_down = _mm_loadu_ps((float const*) &r1[idx]);
	    		__m128 vec_ft = _mm_set_ps(functionTable[idxWhole + 6], functionTable[idxWhole + 4],
	    				functionTable[idxWhole + 2], functionTable[idxWhole]);

	    		__m128 vec_s1 = _mm_add_ps(vec_left, vec_up);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_right);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_down);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_ft);
	    		vec_s1 = _mm_mul_ps(vec_s1, vec_0_25);
	    		_mm_storeu_ps(&s1[idx], vec_s1);

	    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_s1, vec_old));
	    		diff += sse_sum(vec_diff);
	    	}

	    	for (; i < size - 1; i += 2) { // do the remainder sequentially
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		float old_val = s1[idx];

	    		s1[idx] = r1[idx - halfSize] // links
		                  + r1[idx - 1] // oben
		                  + r1[idx + halfSize] // rechts
		                  + r1[idx] // unten
		                  + functionTable[idxWhole];
		        s1[idx] *= 0.25;

		        diff += fabsf(s1[idx] - old_val);
	    	}

	    	for (i = 1; i < size - 8; i += 8) {
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
	    		__m128 vec_old = _mm_loadu_ps((float const*) &s1[idx]);

	    		__m128 vec_left = _mm_loadu_ps((float const*) &r1[idx - halfSize]);
	    		__m128 vec_up = _mm_loadu_ps((float const*) &r1[idx]);
	    		__m128 vec_right = _mm_loadu_ps((float const*) &r1[idx + halfSize]);
	    		__m128 vec_down = _mm_loadu_ps((float const*) &r1[idx + 1]);
	    		__m128 vec_ft = _mm_set_ps(functionTable[idxWhole + 6], functionTable[idxWhole + 4],
	    				functionTable[idxWhole + 2], functionTable[idxWhole]);

	    		__m128 vec_s1 = _mm_add_ps(vec_left, vec_up);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_right);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_down);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_ft);
	    		vec_s1 = _mm_mul_ps(vec_s1, vec_0_25);
	    		_mm_storeu_ps(&s1[idx], vec_s1);

	    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_s1, vec_old));
	    		diff += sse_sum(vec_diff);
	    	}

	    	for (; i < size - 2; i += 2) { // do the remainder sequentially
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
	    		float old_val = s1[idx];

			s1[idx] = r1[idx - halfSize] // links
			          + r1[idx] // oben
			          + r1[idx + halfSize] // rechts
			          + r1[idx + 1] // unten
			          + functionTable[idxWhole];
			s1[idx] *= 0.25;
			diff += fabsf(s1[idx] - old_val);
	    	}
        }

        if (diff / (size * size) < TOL) break;
    }

    #pragma omp parallel for private(i, j)
    for (j = 0; j < size-1; j+=2)
    {
        for (i = 0; i < halfSize; ++i) // for even j
        {
            int baseIdx = j * halfSize;
            int idx = baseIdx + i;
            int idxRot = 2 * idx;
            int idxSchwarz = 2 * idx + 1;
            gaussSeidelResult[idxRot] = r1[idx];
            gaussSeidelResult[idxSchwarz] = s1[idx];
        }
        for (i = 0; i < halfSize; ++i) // for odd j+1
        {
            int baseIdx = (j+1) * halfSize;
            int idx = baseIdx + i;
            int idxRot = 2 * idx + 1;
            int idxSchwarz = 2 * idx;
            gaussSeidelResult[idxRot] = r1[idx];
            gaussSeidelResult[idxSchwarz] = s1[idx];
        }
    }

    free(r1);
    free(s1);
}

// For size % 2 == 1
void gaussSeidelRotSchwarzOdd(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    int i1, i, j, k;
    int halfSize = size / 2;
    int halfSizePlus1 = halfSize + 1;

    unsigned int numRedElements = halfSize * halfSize + halfSizePlus1 * halfSizePlus1;
    unsigned int numBlackElements = 2 * halfSize * halfSizePlus1;

    float* r1 = (float*) malloc(numRedElements * sizeof(float));
    float* s1 = (float*) malloc(numBlackElements * sizeof(float));

    // fill the arrays
    #pragma omp parallel for private(j)
    for (j = 0; j < (size * size) - 1; j+=2)
    {
        // even j:
        r1[j/2] = startVector[j];

        // odd j+1:
        s1[j/2] = startVector[j+1];
    }

    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;
        // rote Punkte
        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 1; i < size - 1; i += 2) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		float old_val = r1[idx];

	    		r1[idx] = s1[idx - halfSizePlus1] // links
			        + s1[idx - 1] // oben
			        + s1[idx + halfSize] // rechts
			        + s1[idx] // unten
			        + functionTable[idx * 2];
			r1[idx] *= 0.25;
			diff += fabsf(r1[idx] - old_val);
	    	}

	    	if (j+1 < size - 1) {
		    	for (i = 2; i < size - 2; i += 2) {
		    		const int idxWhole = CO(i,j+1);
		    		const int idx = idxWhole / 2;
		    		float old_val = r1[idx];

		    		r1[idx] = s1[idx - halfSizePlus1] // links
					+ s1[idx - 1] // oben
					+ s1[idx + halfSize] // rechts
					+ s1[idx] // unten
					+ functionTable[idx * 2];
				r1[idx] *= 0.25;
				diff += fabsf(r1[idx] - old_val);
		    	}
	    	}
        }

        // schwarze Punkte
        #pragma omp parallel for private(j,i) reduction(+:diff)
        for (j = 1; j < size - 1; j += 2) {
        	for (i = 2; i < size - 2; i += 2) {
        		const int idxWhole = CO(i,j);
        		const int idx = idxWhole / 2;
        		float old_val = s1[idx];

        		s1[idx] = r1[idx - halfSize] // links
				        + r1[idx] // oben
				        + r1[idx + halfSizePlus1] // rechts
				        + r1[idx + 1] // unten
				        + functionTable[idxWhole];
				s1[idx] *= 0.25;
			diff += fabsf(s1[idx] - old_val);
        	}

        	if (j+1 < size - 1) {
			for (i = 1; i < size - 1; i += 2) {
				const int idxWhole = CO(i,j+1);
				const int idx = idxWhole / 2;
				float old_val = s1[idx];

				s1[idx] = r1[idx - halfSize] // links
						+ r1[idx] // oben
						+ r1[idx + halfSizePlus1] // rechts
						+ r1[idx + 1] // unten
						+ functionTable[idxWhole];
					s1[idx] *= 0.25;
				diff += fabsf(s1[idx] - old_val);
			}
        	}
        }

        if (diff / (size * size) < TOL) break;
    }

    #pragma omp parallel for private(j)
    for (j = 0; j < (size * size) - 1; j+=2)
    {
        // even j:
        gaussSeidelResult[j] = r1[j/2];
        // odd j+1:
        gaussSeidelResult[j+1] = s1[j/2];
    }

    free(r1);
    free(s1);
}

// For size % 2 == 1
void gaussSeidelRotSchwarzOddSSE(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    const __m128 vec_0_25 = _mm_set_ps1(0.25);
    int i1, i, j, k;
    int halfSize = size / 2;
    int halfSizePlus1 = halfSize + 1;

    unsigned int numRedElements = halfSize * halfSize + halfSizePlus1 * halfSizePlus1;
    unsigned int numBlackElements = 2 * halfSize * halfSizePlus1;

    float* r1 = (float*) malloc(numRedElements * sizeof(float));
    float* s1 = (float*) malloc(numBlackElements * sizeof(float));

    // fill the arrays
    #pragma omp parallel for private(j)
    for (j = 0; j < (size * size) - 1; j+=2)
    {
        // even j:
        r1[j/2] = startVector[j];

        // odd j+1:
        s1[j/2] = startVector[j+1];
    }

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;
        // rote Punkte
        #pragma omp parallel for private(j, i) reduction(+:diff)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 1; i < size - 7; i += 8) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		__m128 vec_old = _mm_loadu_ps((float const*) &r1[idx]);

	    		__m128 vec_left = _mm_loadu_ps((float const*) &s1[idx - halfSizePlus1]);
	    		__m128 vec_up = _mm_loadu_ps((float const*) &s1[idx - 1]);
	    		__m128 vec_right = _mm_loadu_ps((float const*) &s1[idx + halfSize]);
	    		__m128 vec_down = _mm_loadu_ps((float const*) &s1[idx]);
	    		__m128 vec_ft = _mm_set_ps(functionTable[idx * 2 + 6], functionTable[idx * 2 + 4],
	    				functionTable[idx * 2 + 2], functionTable[idx * 2]);

	    		__m128 vec_r1 = _mm_add_ps(vec_left, vec_up);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_right);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_down);
	    		vec_r1 = _mm_add_ps(vec_r1, vec_ft);
	    		vec_r1 = _mm_mul_ps(vec_r1, vec_0_25);
	    		_mm_storeu_ps(&r1[idx], vec_r1);

	    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_r1, vec_old));
	    		diff += sse_sum(vec_diff);
	    	}

	    	for (; i < size - 1; i += 2) { // do the remainder sequentially
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		float old_val = r1[idx];

	    		r1[idx] = s1[idx - halfSizePlus1] // links
			        + s1[idx - 1] // oben
			        + s1[idx + halfSize] // rechts
			        + s1[idx] // unten
			        + functionTable[idx * 2];
			r1[idx] *= 0.25;
			diff += fabsf(r1[idx] - old_val);
	    	}

	    	if (j+1 < size - 1) {
		    	for (i = 2; i < size - 8; i += 8) { // TODO: This is massive code duplication
				const int idxWhole = CO(i, j + 1);
		    		const int idx = idxWhole / 2;
		    		__m128 vec_old = _mm_loadu_ps((float const*) &r1[idx]);

		    		__m128 vec_left = _mm_loadu_ps((float const*) &s1[idx - halfSizePlus1]);
		    		__m128 vec_up = _mm_loadu_ps((float const*) &s1[idx - 1]);
		    		__m128 vec_right = _mm_loadu_ps((float const*) &s1[idx + halfSize]);
		    		__m128 vec_down = _mm_loadu_ps((float const*) &s1[idx]);
		    		__m128 vec_ft = _mm_set_ps(functionTable[idx * 2 + 6], functionTable[idx * 2 + 4],
		    				functionTable[idx * 2 + 2], functionTable[idx * 2]);

		    		__m128 vec_r1 = _mm_add_ps(vec_left, vec_up);
		    		vec_r1 = _mm_add_ps(vec_r1, vec_right);
		    		vec_r1 = _mm_add_ps(vec_r1, vec_down);
		    		vec_r1 = _mm_add_ps(vec_r1, vec_ft);
		    		vec_r1 = _mm_mul_ps(vec_r1, vec_0_25);
		    		_mm_storeu_ps(&r1[idx], vec_r1);

		    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_r1, vec_old));
	    			diff += sse_sum(vec_diff);
		    	}
		    	for (; i < size - 2; i += 2) { // do the remainder sequentially
		    		const int idxWhole = CO(i,j+1);
		    		const int idx = idxWhole / 2;
		    		float old_val = r1[idx];

		    		r1[idx] = s1[idx - halfSizePlus1] // links
					+ s1[idx - 1] // oben
					+ s1[idx + halfSize] // rechts
					+ s1[idx] // unten
					+ functionTable[idx * 2];
				r1[idx] *= 0.25;
				diff += fabsf(r1[idx] - old_val);
		    	}
	    	}
        }

        // schwarze Punkte
        #pragma omp parallel for private(j,i) reduction(+:diff)
        for (j = 1; j < size - 1; j += 2) {
        	for (i = 2; i < size - 8; i += 8) {
        		const int idxWhole = CO(i,j);
        		const int idx = idxWhole / 2;
        		__m128 vec_old = _mm_loadu_ps((float const*) &s1[idx]);

        		__m128 vec_left = _mm_loadu_ps((float const*) &r1[idx - halfSize]);
	    		__m128 vec_up = _mm_loadu_ps((float const*) &r1[idx]);
	    		__m128 vec_right = _mm_loadu_ps((float const*) &r1[idx + halfSizePlus1]);
	    		__m128 vec_down = _mm_loadu_ps((float const*) &r1[idx + 1]);
	    		__m128 vec_ft = _mm_set_ps(functionTable[idxWhole + 6], functionTable[idxWhole + 4],
	    				functionTable[idxWhole + 2], functionTable[idxWhole]);

	    		__m128 vec_s1 = _mm_add_ps(vec_left, vec_up);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_right);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_down);
	    		vec_s1 = _mm_add_ps(vec_s1, vec_ft);
	    		vec_s1 = _mm_mul_ps(vec_s1, vec_0_25);
	    		_mm_storeu_ps(&s1[idx], vec_s1);

	    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_s1, vec_old));
	    		diff += sse_sum(vec_diff);
        	}

        	for (; i < size - 2; i += 2) { // do the remainder sequentially
        		const int idxWhole = CO(i,j);
        		const int idx = idxWhole / 2;
        		float old_val = s1[idx];

        		s1[idx] = r1[idx - halfSize] // links
				        + r1[idx] // oben
				        + r1[idx + halfSizePlus1] // rechts
				        + r1[idx + 1] // unten
				        + functionTable[idxWhole];
				s1[idx] *= 0.25;
			diff += fabsf(s1[idx] - old_val);
        	}

        	if (j+1 < size - 1) {
			for (i = 1; i < size - 7; i += 8) { // TODO: This is massive code duplication
				const int idxWhole = CO(i,j+1);
				const int idx = idxWhole / 2;
				__m128 vec_old = _mm_loadu_ps((float const*) &s1[idx]);

				__m128 vec_left = _mm_loadu_ps((float const*) &r1[idx - halfSize]);
		    		__m128 vec_up = _mm_loadu_ps((float const*) &r1[idx]);
		    		__m128 vec_right = _mm_loadu_ps((float const*) &r1[idx + halfSizePlus1]);
		    		__m128 vec_down = _mm_loadu_ps((float const*) &r1[idx + 1]);
		    		__m128 vec_ft = _mm_set_ps(functionTable[idxWhole + 6], functionTable[idxWhole + 4],
		    				functionTable[idxWhole + 2], functionTable[idxWhole]);

		    		__m128 vec_s1 = _mm_add_ps(vec_left, vec_up);
		    		vec_s1 = _mm_add_ps(vec_s1, vec_right);
		    		vec_s1 = _mm_add_ps(vec_s1, vec_down);
		    		vec_s1 = _mm_add_ps(vec_s1, vec_ft);
		    		vec_s1 = _mm_mul_ps(vec_s1, vec_0_25);
		    		_mm_storeu_ps(&s1[idx], vec_s1);

		    		__m128 vec_diff = abs_ps(_mm_sub_ps(vec_s1, vec_old));
	    			diff += sse_sum(vec_diff);
			}

			for (; i < size - 1; i += 2) { // do the remainder sequentially
				const int idxWhole = CO(i,j+1);
				const int idx = idxWhole / 2;
				const float old_val = s1[idx];
				s1[idx] = r1[idx - halfSize] // links
						+ r1[idx] // oben
						+ r1[idx + halfSizePlus1] // rechts
						+ r1[idx + 1] // unten
						+ functionTable[idxWhole];
					s1[idx] *= 0.25;
				diff += fabsf(s1[idx] - old_val);
			}
        	}
        }

        if (diff / (size * size) < TOL) break;
    }

    #pragma omp parallel for private(j)
    for (j = 0; j < (size * size) - 1; j+=2)
    {
        // even j:
        gaussSeidelResult[j] = r1[j/2];
        // odd j+1:
        gaussSeidelResult[j+1] = s1[j/2];
    }

    free(r1);
    free(s1);
}

bool compare(float* m1,float* m2)
{
    bool equals = true;
    int i =0;
    int j=0;

    for (i = 0; i < (size*size); i++)
    {
        if (fabsf(m1[i]-m2[i]) >= EPSILON)
        {

            equals=false;
        }
        //  printf(" %f ", m1[i*size+j]-m2[i*size+j]);

        //printf("\n");
    }
    return equals;
}
void computeFunctionTable(float h, float* functionTable)
{
    int i, j;
    float hSquared = h * h;
    #pragma omp parallel for private(i, j)
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            float x1 = i * h;
            float y1 = j * h;
            functionTable[CO(i, j)] = hSquared * f1(x1, y1);
        }
    }
}

float f1(float x, float y)
{
    return (32 * (x * (1 - x) + y * (1 - y)));
}

void printResultMatrix(const float* matrix)
{
    int i, j;
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            float val = matrix[CO(i, j)];
            if (j == 0)
            {
                printf("%.3f", val);
            }
            else if (j < size - 1)
            {
                printf(" %.3f ", val);
            }
            else
            {
                printf(" %.3f\n", val);
            }
        }
    }
}

void printAnalyticalResult(float h)
{
    int i, j;
    for (i = 0; i < size; ++i)
    {
        float x = i * h;
        for (j = 0; j < size; ++j)
        {
            float y = j * h;
            float val;
            if (i == 0 || j == 0 || i == size - 1 || j == size - 1)
            {
                val = 0; // Randbedingung
            }
            else
            {
                val = 16 * x * (1 - x) * y * (1 - y);
            }

            if (j == 0)
            {
                printf("%.3f", val);
            }
            else if (j < size - 1)
            {
                printf(" %.3f ", val);
            }
            else
            {
                printf(" %.3f\n", val);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: ./parallelDG SIZE\n (SIZE stands for the number of grid points in each dimension)\n");
        return -1;
    }

    size = atoi(argv[1]);
    if (size <= 1)
    {
        printf("Error! Size must be greater than 1.\n");
        return -1;
    }
    float h = 1 / (float) (size - 1);
    // Precompute function f
    float* precomputedF = malloc(size * size * sizeof(float));
    computeFunctionTable(h, precomputedF);
    // Analytical result
    float* analyticalResult = (float*) malloc(size * size * sizeof(float));
    // Create random start vector with zeroes at the proper positions (at the border)
    float* startVector = (float*) malloc(size * size * sizeof(float));
    srand(time(NULL));
    int i, j;

    #pragma omp parallel for private(i, j)
    for (i = 0; i < size; ++i)
    {
    	float x = i * h;
        for (j = 0; j < size; ++j)
        {
            float y = j * h;
            float val;
            if (i == 0 || j == 0 || i == size - 1 || j == size - 1)
            {
                val = 0; // Randbedingung
            }
            else
            {
                val = 16 * x * (1 - x) * y * (1 - y);
            }
            analyticalResult[CO(i,j)] = val;

            if (i != 0 && j != 0 && i != size - 1 && j != size - 1)
            {
                val = (float) rand() / RAND_MAX;
            }
            startVector[CO(i, j)] = val;
        }
    }

    double start, end;

    int repeats = 1;

    // Call Jacobi Sequential
    float* jacobiSequentialResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	jacobiSerial(startVector, h, precomputedF, jacobiSequentialResult);
    }
    end = get_wall_time();
    printf("Execution time Jacobi Sequential: %.3f seconds\n", (end - start) / repeats);
    bool correct=true;
    correct=compare(jacobiSequentialResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");


    // Call Jacobi
    float* jacobiResult = malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	jacobi(startVector, h, precomputedF, jacobiResult);
    }
    end = omp_get_wtime();
    printf("Execution time Jacobi: %.3f seconds\n", (end - start) / repeats);
    correct=compare(jacobiResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");

    // Call Jacobi SSE
    float* jacobiSSEResult = malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	jacobiSSE(startVector, h, precomputedF, jacobiSSEResult);
    }
    end = omp_get_wtime();
    printf("Execution time Jacobi SSE: %.3f seconds\n", (end - start) / repeats);
    correct=compare(jacobiSSEResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");

    // Call Gauss-Seidel
    float* gaussSeidelResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidel(startVector, h, precomputedF, gaussSeidelResult);
    }
    end = get_wall_time();
    printf("Execution time Gauss-Seidel: %.3f seconds\n", (end - start) / repeats);
    correct=compare(gaussSeidelResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");

    // Call Gauss-Seidel Naiv
    float* gaussSeidelNaivResult = malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelNaiv(startVector, h, precomputedF, gaussSeidelNaivResult);
    }
    end = omp_get_wtime();
    printf("Execution time Gauss-Seidel Naiv: %.3f seconds\n", (end - start) / repeats);
    correct=true;
    correct=compare(gaussSeidelNaivResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");

    // Call Gauss-Seidel Rot-Schwarz
    float* gaussSeidelRotSchwarzResult = malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelRotSchwarz(startVector, h, precomputedF, gaussSeidelRotSchwarzResult);
    }
    end = omp_get_wtime();
    printf("Execution time Gauss-Seidel Rot-Schwarz: %.3f seconds\n", (end - start) / repeats);
    correct=compare(gaussSeidelRotSchwarzResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");

    // Call Gauss-Seidel Rot-Schwarz SSE
    float* gaussSeidelRotSchwarzSSEResult = malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelRotSchwarzSSE(startVector, h, precomputedF, gaussSeidelRotSchwarzSSEResult);
    }
    end = omp_get_wtime();
    printf("Execution time Gauss-Seidel Rot-Schwarz SSE: %.3f seconds\n", (end - start) / repeats);
    correct=compare(gaussSeidelRotSchwarzSSEResult, analyticalResult);
    printf("  is it correct: %s  \n" ,(correct)?"true":"false");

    //Call Gaus Seidel Wavefront
    float* gaussSeidelWavefrontResult= malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelWavefront(startVector, h, precomputedF, gaussSeidelWavefrontResult);
    }
    end = omp_get_wtime();
    printf("Execution time Gauss-Seidel Wavefront: %.3f seconds\n", (end - start) / repeats);
    correct=compare(gaussSeidelWavefrontResult, analyticalResult);
    printf("  is it correct: %s \n" ,(correct)?"true":"false");

    //Call Gaus Seidel Wavefront Cache
    float* gaussSeidelWavefrontCacheResult= malloc(size * size * sizeof(float));
    start = omp_get_wtime();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelWavefrontCache(startVector, h, precomputedF, gaussSeidelWavefrontCacheResult);
    }
    end = omp_get_wtime();
    printf("Execution time Gauss-Seidel WavefrontCache: %.3f seconds\n", (end - start) / repeats);
    correct=compare(gaussSeidelWavefrontCacheResult, analyticalResult);
    printf("  is it correct: %s \n" ,(correct)?"true":"false");

    // TODO: The following is just debug code. Remove afterwards.
    /*printf("\nFunctionTable:\n");
    printResultMatrix(precomputedF);
    printf("\nStartvektor:\n");
    printResultMatrix(startVector);
    printf("\nErgebnis Jacobi-Verfahren:\n");
    printResultMatrix(jacobiResult);
    printf("\nErgebnis Jacobi-Verfahren SSE:\n");
    printResultMatrix(jacobiSSEResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren:\n");
    printResultMatrix(gaussSeidelResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Rot-Schwarz:\n");
    printResultMatrix(gaussSeidelRotSchwarzResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Rot-Schwarz SSE:\n");
    printResultMatrix(gaussSeidelRotSchwarzSSEResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Wavefront:\n");
    printResultMatrix(gaussSeidelWavefrontResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Wavefront Cache:\n");
    printResultMatrix(gaussSeidelWavefrontCacheResult);
    printf("\nErgebnis analytisch:\n");
    printAnalyticalResult(h);*/

    free(gaussSeidelWavefrontResult);
    free(gaussSeidelWavefrontCacheResult);
    free(jacobiResult);
    free(jacobiSSEResult);
    free(gaussSeidelResult);
    free(gaussSeidelRotSchwarzResult);
    free(gaussSeidelRotSchwarzSSEResult);
    free(startVector);
    free(precomputedF);
    free(analyticalResult);
    return 0;
}

void gaussSeidelWavefront(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{

    float* a0 = (float*) malloc(size * size * sizeof(float));
    float* a1 = (float*) malloc(size * size * sizeof(float));


    int i;
    #pragma omp parallel for
    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    int k = 0;

    for (k = 0; k < MAX_ITERATIONS; k++)
    {
	float diff = 0;
        float* temp = a0;
        a0 = a1;
        a1 = temp;

        int currentEle = 0;
        int border = 0;
        int durchlauf;
        for (durchlauf = 0; durchlauf<size+size-1-4 ; durchlauf++) //-1 weil diagonalen zahl size+size-1, -4 weil 4 diagonalen wegfallen
        {
            if (durchlauf > (size -1-2))//-1 weil fängt bei 0 an -2 weil ersten 2 diagonalen rand sind
            {
                currentEle--;
                border++;
            }
            else
            {
                currentEle++;
            }
            int i = 0;
           #pragma omp parallel for firstprivate(durchlauf,border,currentEle,k) reduction(+:diff)
            for (i = 0; i < currentEle; i++)
            {

                //  printf("(%i %i %i %i)",index-1-newsize,index-newsize,index+1+newsize,index+newsize);
                int index= (durchlauf - border - i+1)* size +( i + border+1);
                a1[index] = 0.25   * (a1[index-1]
                                      + a1[index-size]
                                      + a0[index+1]
                                      + a0[index+size]
                                      +  functionTable[index]);

                diff += fabsf(a1[index] - a0[index]);
      /*          printf("(%i)",index);
                printf("(%i %i %i %i)",index-1,index-size,index+1,index+size);
                printf("(%f %f %f %f %f %f) \n",a1[index], a1[index-1],a1[index-size], a0[index+1], a0[index+size]    ,  functionTable[index]);*/
                /*    a1[(durchlauf - border - i+1)* size +( i + border+1)] = 0.25 //+1 jeweils für den rand dei anderen indexe sind algorythmus relevant
                            * (a1[(durchlauf - border - i+1)* size  + (i + border - 1+1)]
                               + a1[(durchlauf - border - 1 - i+1)* size  + (i + border+1)]
                               + a0[(durchlauf - border + 1 - i+1)* size  + (i + border+1)]
                               + a0[(durchlauf - border - i+1)* size  + (i + border + 1+1)]
                               +  functionTable[(durchlauf - border - i+1) * size + (i + border+1)]); */
            }
        }

	if (diff / (size * size) < TOL) break;
    }
    #pragma omp parallel for
    for (i = 0; i < size * size; i++)
    {
        gaussSeidelResult[i] = a1[i];

    }
    free(a0);
    free(a1);


}
void gaussSeidelWavefrontCache(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    int newsize= (2*size-1);
    float* a0 = (float*) malloc((2*size-1) * newsize  * sizeof(float));
    float* a1 = (float*) malloc((2*size-1) * newsize  * sizeof(float));

    int i;
    int currentEle = 0;
    int border = 0;
    int durchlauf;
    //kopieren
      #pragma omp parallel for firstprivate(border,currentEle)private(i,durchlauf)
    for (durchlauf = 0; durchlauf<newsize ; durchlauf++) //-1 weil diagonalen zahl size+size-1
    {

        if (durchlauf > (size -1))//-1 weil fängt bei 0 an -2 weil ersten 2 diagonalen rand sind
        {
            currentEle--;
            border++;
        }
        else
        {
            currentEle++;
        }
        int i = 0;

        for (i = 0; i < currentEle; i++)
        {
            //(durchlauf - border - i+1)* size +( i + border+1);
            int indexZu=durchlauf*newsize+i;
            int indexVon=(durchlauf - border - i)* (size) +( i + border);
            // printf("%i %i %i \n" ,durchlauf,border,i);
            a0[indexZu]= startVector[indexVon];
            a1[indexZu]= startVector[indexVon];

            //printf("%f ",a1[durchlauf*size+i]);
         //  printf("%i<-%i ",indexZu,indexVon);
            // printf(" %i ",indexZu);
            /* printf("%i ",i);
             printf("%i ",durchlauf); */
        }
       // printf("\n");
    }





//arbeiten


//todo abbruchbedingung Max itera wieder ienführen
    int k = 0;

    for (k = 0; k < MAX_ITERATIONS; k++)
    {
        float diff=0;
        float* temp = a0;
        a0 = a1;
        a1 = temp;

        currentEle = 0;
        border = 0;
      //  printf("arbeiten..:\n");

        for (durchlauf = 2; durchlauf<newsize-2 ; durchlauf++) //-1 weil diagonalen zahl size+size-1, -4 weil 4 diagonalen wegfallen
        {


            if (durchlauf > (size -1))//-1 weil fängt bei 0 an -2 weil ersten 2 diagonalen rand sind
            {

                currentEle--;
                border++;
            }

            else
            {
                currentEle++;
            }
            int i = 0;
                 #pragma omp parallel for firstprivate(durchlauf,border,currentEle,k)reduction(+:diff)
            for (i = 1; i < currentEle+1; i++)
            {
                /* int indexZu=durchlauf*newsize+i;
            int indexVon=(durchlauf - border - i)* (size) +( i + border);
                int indexZu=durchlauf*(size+size-1)+i;
                int indexVon=(durchlauf - border - i)* (size) +( i + border);
                   int index= (durchlauf - border - i+1)* size +( i + border+1);
                a1[index] = 0.25   * (a1[index-1]
                                      + a1[index-size]
                                      + a0[index+1]
                                      + a0[index+size]
                                      +  functionTable[index]);
                */
           /*     int indexVon=(durchlauf - border - i+1)* (size) +( i + border+1);
                int index= (durchlauf+2 )* newsize +( i +1);//+2 wegen erste 2 diagonalen müll +1 weil die zahlen 1 drinne stehen
                int a11=index-1-newsize+switchIt;
                int a12 = index-newsize+switchIt;
                int a01=index+newsize-border;
                int a02=index+1+newsize-border-switchA1; */
             int hack2=0;
        int hack=0;
                int indexZu=durchlauf*newsize+i;
                 hack=(int)(durchlauf/(size-1));
                 hack2=(int)(durchlauf/(size));

                int index=(durchlauf - border - i+1)* size +( i + border+1);

                int a11=(durchlauf-1)*newsize+i-1+hack2;//+(border);//;+hack2;
                int a12=(durchlauf-1)*newsize+i+hack2;// +(border);//;+hack2;

                int a01=(durchlauf+1)*newsize+i+1-hack;
                int a02=(durchlauf+1)*newsize+i-hack;
                int indexVon=(durchlauf - border - i)* (size) +( i + border);

                a1[indexZu] = 0.25 * (a1[a11]
                                    + a1[a12]
                                    + a0[a01]
                                    + a0[a02]
                                    +  functionTable[indexVon]);
                                     diff += fabsf(a1[indexZu] - a0[indexZu]);
                //  printf("%i ",index);
                //  printf("(%i asdf %i %i %i)",indexVon,durchlauf,border,i);
         //       printf("(a1 %i F %i Index %i)",indexZu,indexVon,index);
           //     printf("(durchlauf %i border  %i hack %i cuele %i size %i)",durchlauf,border,hack,currentEle,size);
             //   printf("(%i %i %i %i) \n",a11,a12,a01,a02);
              //  printf("(%f %f %f %f %f %f) \n",a1[indexVon],a1[a11],a1[a12],a0[a02],a0[a01],functionTable[indexVon]);

            }


        }
if (diff / (size * size) < TOL) break;

    }
    //zurückopieren

    i=0;
    currentEle = 0;
    border = 0;
 //   printf("zurückopieren..:\n");
    //kopieren
    #pragma omp parallel for firstprivate(border,currentEle) private(i,durchlauf)
    for (durchlauf = 0; durchlauf<newsize ; durchlauf++) //-1 weil diagonalen zahl size+size-1
    {

        if (durchlauf > (size -1))//-1 weil fängt bei 0 an -2 weil ersten 2 diagonalen rand sind
        {
            currentEle--;
            border++;
        }
        else
        {
            currentEle++;
        }
        int i = 0;
        //     #pragma omp parallel for firstprivate(durchlauf,border,currentEle,k)
        for (i = 0; i < currentEle; i++)
        {
            int indexVon=durchlauf*newsize+i;
            int indexZu=(durchlauf - border - i)* (size) +( i + border);
            gaussSeidelResult[indexZu]= a1[indexVon];
      //      printf("(%i<- %i)",indexZu,indexVon);

        }
     //   printf("\n");
    }



    free(a0);
    free(a1);


}


void gaussSeidelNaiv(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult)
{
    int i, j, k;

    float* array0 = (float*) malloc(size * size * sizeof(float));
    float* array1 = (float*) malloc(size * size * sizeof(float));

    float* a0 = array0; // last iteration
    float* a1 = array1; // current iteration
    #pragma omp parallel for
    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
    	float diff = 0;
        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;
        #pragma omp parallel for private(j, i) collapse(2) reduction(+:diff)
        for (j = 1; j < size - 1; j++)
        {
            for (i = 1; i < size - 1; i++)
            {
                a1[CO(i,j)] = a1[CO(i, j - 1)]
                              + a1[CO(i - 1, j)]
                              + a0[CO(i, j + 1)]
                              + a0[CO(i + 1, j)]
                              + functionTable[CO(i, j)];
                a1[CO(i,j)] *= 0.25;

                diff += fabsf(a1[CO(i,j)] - a0[CO(i,j)]);
            }
        }
        if (diff / (size * size) < TOL) break;
    }
    #pragma omp parallel for
    for (i = 0; i < size * size; ++i)
    {
        gaussSeidelResult[i] = a1[i];
    }
    free(a0);
    free(a1);
}

