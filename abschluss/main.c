#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

/*
Compile with: gcc -O2 -fopenmp -march=native main.c -o main
*/

float f1(float x, float y);
void jacobiSerial(const float* startVector, float h, const float* functionTable, float* jacobiResult);
void jacobi(const float* startVector, float h, const float* functionTable, float* jacobiResult);
void gaussSeidel(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzEven(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarzOdd(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelNaiv(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarz(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void computeFunctionTable(float h, float* functionTable);
void printResultMatrix(const float* matrix);
void printAnalyticalResult(float h);
void gaussSeidelWavefront(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelWavefrontCache(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
bool compare(float* m1,float* m2);

int size;

#define CO(i,j) ( (j) * (size) + (i) )

const static int MAX_ITERATIONS = 10000;
const static float EPSILON = 0.00000001;

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

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
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
            }
        }
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

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;
        #pragma omp parallel for private(j, i) collapse(2)
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
            }
        }
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

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;
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


            }
        }
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

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {   
        // rote Punkte
        #pragma omp parallel for private(j, i)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 1; i < size - 2; i += 2) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
			r1[idx] = s1[idx - halfSize] // links
		                  + s1[idx] // oben
		                  + s1[idx + halfSize] // rechts
		                  + s1[idx + 1] // unten
		                  + functionTable[idxWhole];
		        r1[idx] *= 0.25;
	    	}
	    	
	    	for (i = 2; i < size - 1; i += 2) {
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
			
			r1[idx] = s1[idx - halfSize] // links
			          + s1[idx - 1] // oben
			          + s1[idx + halfSize] // rechts
			          + s1[idx] // unten
			          + functionTable[idxWhole];
			r1[idx] *= 0.25;
	    	}
        }
        
        // schwarze Punkte
        #pragma omp parallel for private(j, i)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 2; i < size - 1; i += 2) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		s1[idx] = r1[idx - halfSize] // links
		                  + r1[idx - 1] // oben
		                  + r1[idx + halfSize] // rechts
		                  + r1[idx] // unten
		                  + functionTable[idxWhole];
		        s1[idx] *= 0.25;
	    	}
	    	
	    	for (i = 1; i < size - 2; i += 2) {
	    		const int idxWhole = CO(i,j+1);
	    		const int idx = idxWhole / 2;
			
			s1[idx] = r1[idx - halfSize] // links
			          + r1[idx] // oben
			          + r1[idx + halfSize] // rechts
			          + r1[idx + 1] // unten
			          + functionTable[idxWhole];
			s1[idx] *= 0.25;
	    	}
        }
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

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {   
        // rote Punkte, neuer Versuch
        #pragma omp parallel for private(j, i)
        for (j = 1; j < size - 1; j+=2)
        {
	    	for (i = 1; i < size - 1; i += 2) {
	    		const int idxWhole = CO(i,j);
	    		const int idx = idxWhole / 2;
	    		r1[idx] = s1[idx - halfSizePlus1] // links
			        + s1[idx - 1] // oben
			        + s1[idx + halfSize] // rechts
			        + s1[idx] // unten
			        + functionTable[idx * 2];
			r1[idx] *= 0.25;
	    	}
	    	
	    	if (j+1 < size - 1) {
		    	for (i = 2; i < size - 2; i += 2) {
		    		const int idxWhole = CO(i,j+1);
		    		const int idx = idxWhole / 2;
		    		r1[idx] = s1[idx - halfSizePlus1] // links
					+ s1[idx - 1] // oben
					+ s1[idx + halfSize] // rechts
					+ s1[idx] // unten
					+ functionTable[idx * 2];
				r1[idx] *= 0.25;
		    	}
	    	}
        }
        
        // schwarze Punkte, und wieder neuer Versuch
        #pragma omp parallel for private(j,i)
        for (j = 1; j < size - 1; j += 2) {
        	for (i = 2; i < size - 2; i += 2) {
        		const int idxWhole = CO(i,j);
        		const int idx = idxWhole / 2;
        		s1[idx] = r1[idx - halfSize] // links
				        + r1[idx] // oben
				        + r1[idx + halfSizePlus1] // rechts
				        + r1[idx + 1] // unten
				        + functionTable[idxWhole];
				s1[idx] *= 0.25;
        	}
        	
        	if (j+1 < size - 1) {
			for (i = 1; i < size - 1; i += 2) {
				const int idxWhole = CO(i,j+1);
				const int idx = idxWhole / 2;
				s1[idx] = r1[idx - halfSize] // links
						+ r1[idx] // oben
						+ r1[idx + halfSizePlus1] // rechts
						+ r1[idx + 1] // unten
						+ functionTable[idxWhole];
					s1[idx] *= 0.25;
			}
        	}
        }
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

    for(i=0; i<(size*size); i++)
    {
        if((m1[i]-m2[i]) * (m1[i]-m2[i]) >= EPSILON)
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
    // Create random start vector with zeroes at the proper positions (at the border)
    float* startVector = (float*) malloc(size * size * sizeof(float));
    srand(time(NULL));
    int i, j;

    #pragma omp parallel for private(i, j)
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            float val = 0;
            if (i != 0 && j != 0 && i != size - 1 && j != size - 1)
            {
                val = (float) rand() / RAND_MAX;
            }
            startVector[CO(i, j)] = val;
        }
    }

    double start, end;
    
    int repeats = 1;

    // Call JacobiSerial
    float* jacobiSerialResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	jacobiSerial(startVector, h, precomputedF, jacobiSerialResult);
    }
    end = get_wall_time();
    printf("Execution time JacobiSerial: %.3f seconds\n", (end - start) / repeats);

    // Call Jacobi
    float* jacobiResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	jacobi(startVector, h, precomputedF, jacobiResult);
    }
    end = get_wall_time();
    printf("Execution time Jacobi: %.3f seconds\n", (end - start) / repeats);

    // Call Gauss-Seidel
    float* gaussSeidelResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidel(startVector, h, precomputedF, gaussSeidelResult);
    }
    end = get_wall_time();
    printf("Execution time Gauss-Seidel: %.3f seconds\n", (end - start) / repeats);

    // Call Gauss-Seidel
    float* gaussSeidelNaivResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelNaiv(startVector, h, precomputedF, gaussSeidelNaivResult);
    }
    end = get_wall_time();
    printf("Execution time Gauss-Seidel Naiv: %.3f seconds", (end - start) / repeats);
    bool correct=true;
    correct=compare(gaussSeidelNaivResult,gaussSeidelResult);
    printf("is it correct: %s  \n" ,(correct)?"true":"false");

    // Call Gauss-Seidel Rot-Schwarz
    float* gaussSeidelRotSchwarzResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelRotSchwarz(startVector, h, precomputedF, gaussSeidelRotSchwarzResult);
    }
    end = get_wall_time();
    printf("Execution time Gauss-Seidel Rot-Schwarz: %.3f seconds", (end - start) / repeats);
    correct=compare(gaussSeidelRotSchwarzResult,gaussSeidelResult);
    printf("is it correct: %s  \n" ,(correct)?"true":"false");

    //Call Gaus Seidel Wavefront
    float* gaussSeidelWavefrontResult= malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelWavefront(startVector, h, precomputedF, gaussSeidelWavefrontResult);
    }
    end = get_wall_time();
    printf("Execution time Gauss-Seidel Wavefront: %.3f seconds ", (end - start) / repeats);
    correct=compare(gaussSeidelWavefrontResult,gaussSeidelResult);
    printf("is it correct: %s \n" ,(correct)?"true":"false");

    //Call Gaus Seidel Wavefront Cache
    float* gaussSeidelWavefrontCacheResult= malloc(size * size * sizeof(float));
    start = get_wall_time();
    for (i = 0; i < repeats; ++i) {
    	gaussSeidelWavefrontCache(startVector, h, precomputedF, gaussSeidelWavefrontCacheResult);
    }
    end = get_wall_time();
    printf("Execution time Gauss-Seidel WavefrontCache: %.3f seconds ", (end - start) / repeats);
    correct=compare(gaussSeidelWavefrontCacheResult,gaussSeidelResult);
    printf("is it correct: %s \n" ,(correct)?"true":"false");

    // TODO: The following is just debug code. Remove afterwards.
    /*printf("\nFunctionTable:\n");
    printResultMatrix(precomputedF);
    printf("\nStartvektor:\n");
    printResultMatrix(startVector);
    printf("\nErgebnis Jacobi-Verfahren:\n");
    printResultMatrix(jacobiResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren:\n");
    printResultMatrix(gaussSeidelResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Rot-Schwarz:\n");
    printResultMatrix(gaussSeidelRotSchwarzResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Wavefront:\n");
    printResultMatrix(gaussSeidelWavefrontResult);
    printf("\nErgebnis Gauss-Seidel-Verfahren Wavefront Cache:\n");
    printResultMatrix(gaussSeidelWavefrontCacheResult);
    printf("\nErgebnis analytisch:\n");
    printAnalyticalResult(h);*/

    free(gaussSeidelWavefrontResult);
    free(gaussSeidelWavefrontCacheResult);
    free(jacobiResult);
    free(gaussSeidelResult);
    free(gaussSeidelRotSchwarzResult);
    free(startVector);
    free(precomputedF);
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

//todo abbruchbedingung Max itera wieder ienführen


    for (k = 0; k < MAX_ITERATIONS; k++)
    {

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
           #pragma omp parallel for firstprivate(durchlauf,border,currentEle,k)
            for (i = 0; i < currentEle; i++)
            {

                //  printf("(%i %i %i %i)",index-1-newsize,index-newsize,index+1+newsize,index+newsize);
                int index= (durchlauf - border - i+1)* size +( i + border+1);
                a1[index] = 0.25   * (a1[index-1]
                                      + a1[index-size]
                                      + a0[index+1]
                                      + a0[index+size]
                                      +  functionTable[index]);
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
      #pragma omp parallel for firstprivate(border,currentEle)private(i)
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
        //    printf("%i<-%i ",indexZu,indexVon);
            // printf(" %i ",indexZu);
            /* printf("%i ",i);
             printf("%i ",durchlauf); */
        }
    //    printf("\n");
    }





//arbeiten


//todo abbruchbedingung Max itera wieder ienführen
    int k = 0;

    for (k = 0; k < MAX_ITERATIONS; k++)
    {

        float* temp = a0;
        a0 = a1;
        a1 = temp;

        currentEle = 0;
        border = 0;
      //  printf("arbeiten..:\n");
        for (durchlauf = 0; durchlauf<newsize-4 ; durchlauf++) //-1 weil diagonalen zahl size+size-1, -4 weil 4 diagonalen wegfallen
        {
            int switchIt=0;
            int switchA1=0;
            if (durchlauf > (size -1-2))//-1 weil fängt bei 0 an -2 weil ersten 2 diagonalen rand sind
            {
            switchIt=1;

                currentEle--;

                border++;
            }
            else if(durchlauf<size-1-2){


                currentEle++;
            }
            else
            {

            switchA1=2;
                currentEle++;
            }
            int i = 0;
                 #pragma omp parallel for firstprivate(durchlauf,border,currentEle,k)
            for (i = 0; i < currentEle; i++)
            {
                /*int indexZu=durchlauf*(size+size-1)+i;
                int indexVon=(durchlauf - border - i)* (size) +( i + border);
                */
                int indexVon=(durchlauf - border - i+1)* (size) +( i + border+1);
                int index= (durchlauf+2 )* newsize +( i +1);//+2 wegen erste 2 diagonalen müll +1 weil die zahlen 1 drinne stehen
                int a11=index-1-newsize+switchIt;
                int a12 = index-newsize+switchIt;
                int a01=index+newsize-border;
                int a02=index+1+newsize-border-switchA1;

                a1[index] = 0.25 * (a1[a11]
                                    + a1[a12]
                                    + a0[a01]
                                    + a0[a02]
                                    +  functionTable[indexVon]);
                //  printf("%i ",index);
                //  printf("(%i asdf %i %i %i)",indexVon,durchlauf,border,i);
          /*      printf("(%i<- %i)",indexVon,index);
                printf("(durchlauf %i border  %i switch %i )\n",durchlauf,border,switchIt);
                printf("(%i %i %i %i)",a11,a12,a01,a02);
                printf("(%f %f %f %f %f %f) \n",a1[index],a1[a11],a1[a12],a0[a02],a0[a01],functionTable[indexVon]); */

            }

        }


    }
    //zurückopieren

    i=0;
    currentEle = 0;
    border = 0;
 //   printf("zurückopieren..:\n");
    //kopieren
    #pragma omp parallel for firstprivate(border,currentEle) private(i)
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
        // swap a0 and a1
        float* temp = a0;
        a0 = a1;
        a1 = temp;
        #pragma omp parallel for private(j, i) collapse(2)
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


            }
        }
    }
    #pragma omp parallel for
    for (i = 0; i < size * size; ++i)
    {
        gaussSeidelResult[i] = a1[i];
    }
    free(a0);
    free(a1);
}

