#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

float f1(float x, float y);
void jacobi(const float* startVector, float h, const float* functionTable, float* jacobiResult);
void gaussSeidel(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void gaussSeidelRotSchwarz(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
void computeFunctionTable(float h, float* functionTable);
void printResultMatrix(const float* matrix);
void printAnalyticalResult(float h);
void gaussSeidelWavefront(const float * startVector, float h, const float* functionTable, float* gaussSeidelResult);
bool compare(float* m1,float* m2);

int size;

#define CO(i,j) ( (j) * (size) + (i) )

static int MAX_ITERATIONS = 10000;

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
    int i1, i, j, k;
    int halfSize = size / 2;
    float* arrayRot0 = (float*) malloc(size * halfSize * sizeof(float));
    float* arrayRot1 = (float*) malloc(size * halfSize * sizeof(float));
    float* arraySchwarz0 = (float*) malloc(size * halfSize * sizeof(float));
    float* arraySchwarz1 = (float*) malloc(size * halfSize * sizeof(float));

    float* r0 = arrayRot0; // last iteration Rot
    float* r1 = arrayRot1; // current iteration Rot
    float* s0 = arraySchwarz0; // last iteration Schwarz
    float* s1 = arraySchwarz1; // current iteration Schwarz

    #pragma omp parallel for private(i, j) collapse(2)
    for (j = 0; j < size; ++j)
    {
        //printf("Spalte %d:\n", j);
        //printf(" IdxRot: ");
        for (i = 0; i < halfSize; ++i)
        {
            int baseIdx = j * halfSize;
            int idx = baseIdx + i;
            int idxRot, idxSchwarz;
            if (j % 2 == 0)
            {
                idxRot = 2 * idx;
                idxSchwarz = 2 * idx + 1;
            }
            else
            {
                idxRot = 2 * idx + 1;
                idxSchwarz = 2 * idx;
            }

            //printf("%d ", idxRot);
            arrayRot0[idx] = startVector[idxRot];
            arrayRot1[idx] = startVector[idxRot];
            arraySchwarz0[idx] = startVector[idxSchwarz];
            arraySchwarz1[idx] = startVector[idxSchwarz];
        }
        //printf("\n");
    }

    /*printf("\narrayRot:\n");
    for (i = 0; i < size * size/2; ++i) {
    	printf("%.3f ", arrayRot0[i]);
    }
    printf("\n");

    printf("\narraySchwarz:\n");
    for (i = 0; i < size * size/2; ++i) {
    	printf("%.3f ", arraySchwarz0[i]);
    }
    printf("\n");*/

    //todo abbruchbedingung
    for (k = 0; k < MAX_ITERATIONS; ++k)
    {
        // swap the pointers for current and last iteration
        float* temp = r0;
        r0 = r1;
        r1 = temp;
        temp = s0;
        s0 = s1;
        s1 = temp;

        // rote Punkte
        #pragma omp parallel for private(i, j) firstprivate(s0, s1) collapse(2)
        for (j = 1; j < size - 1; j++)
        {
            for (i1 = halfSize; i1 < size - 1; i1++)
            {
                int y = j;
                bool offsetRot = (j % 2 == 0);
                const int baseIdx = (j - 1) * halfSize + offsetRot;

                const int idx = i1 + baseIdx;
                r1[idx] = s1[idx - halfSize] // links
                          + s1[idx - offsetRot] // oben
                          + s0[idx + halfSize] // rechts
                          + s0[idx + !offsetRot] // unten
                          + functionTable[idx * 2 + !offsetRot];
                r1[idx] *= 0.25;

                /*printf("  Rot-Index: %d\n", idx);
                printf("    offsetRot: %d\n", offsetRot);
                printf("    Schwarz-Index links: %d enthält: %.3f\n", idx - halfSize, s1[idx - halfSize]);
                printf("    Schwarz-Index oben: %d enthält: %.3f\n", idx - offsetRot, s1[idx - offsetRot]);
                printf("    Schwarz-Index rechts: %d enthält: %.3f\n", idx + halfSize, s0[idx + halfSize]);
                printf("    Schwarz-Index unten: %d enthält: %.3f\n", idx + !offsetRot, s0[idx + !offsetRot]);
                printf("    Zugriff auf f bei: %d\n", idx * 2 + !offsetRot);
                printf("    Ergebnis: %.3f\n", r1[idx]);*/
            }
        }

        // schwarze Punkte
        #pragma omp parallel for private(i, j) firstprivate(r0, r1) collapse(2)
        for (j = 1; j < size - 1; j++)
        {
            for (i1 = halfSize; i1 < size - 1; i1++)
            {
                int y = j;
                bool offsetSchwarz = (j % 2 != 0);
                const int baseIdx = (j - 1) * halfSize + offsetSchwarz;

                const int idx = i1 + baseIdx;
                s1[idx] = r1[idx - halfSize] // links
                          + r1[idx - offsetSchwarz] // oben
                          + r0[idx + halfSize] // rechts
                          + r0[idx + !offsetSchwarz] // unten
                          + functionTable[idx * 2 + !offsetSchwarz];
                s1[idx] *= 0.25;

                /*printf("  Schwarz-Index: %d\n", idx);
                printf("    offsetSchwarz: %d\n", offsetSchwarz);
                printf("    Rot-Index links: %d enthält: %.3f\n", idx - size/2, s1[idx - size/2]);
                printf("    Rot-Index oben: %d enthält: %.3f\n", idx - offsetSchwarz, s1[idx - offsetSchwarz]);
                printf("    Rot-Index rechts: %d enthält: %.3f\n", idx + size/2, s0[idx + size/2]);
                printf("    Rot-Index unten: %d enthält: %.3f\n", idx + !offsetSchwarz, s0[idx + !offsetSchwarz]);
                printf("    Zugriff auf f bei: %d\n", idx * 2 + !offsetSchwarz);
                printf("    Ergebnis: %.3f\n", s1[idx]);*/
            }
        }
    }

    /*printf("\nr1:\n");
    for (i = 0; i < size * size/2; ++i) {
    	printf("%.3f ", r1[i]);
    }
    printf("\n");

    printf("\ns1:\n");
    for (i = 0; i < size * size/2; ++i) {
    	printf("%.3f ", s1[i]);
    }
    printf("\n");*/

    #pragma omp parallel for private(i, j) collapse(2)
    for (j = 0; j < size; ++j)
    {
        for (i = 0; i < halfSize; ++i)
        {
            int baseIdx = j * halfSize;
            int idx = baseIdx + i;
            int idxRot, idxSchwarz;
            if (j % 2 == 0)
            {
                idxRot = 2 * idx;
                idxSchwarz = 2 * idx + 1;
            }
            else
            {
                idxRot = 2 * idx + 1;
                idxSchwarz = 2 * idx;
            }
            gaussSeidelResult[idxRot] = r1[idx];
            gaussSeidelResult[idxSchwarz] = s1[idx];
        }
    }

    free(r0);
    free(r1);
    free(s0);
    free(s1);
}

    bool compare(float* m1,float* m2)
    {
        bool equals = true;
        int i =0;
        printf("\n");
        //#pragma omp parallel for private (i) shared(equals)
        for(i=0; i<(size*size); i++)
        {
       //     printf("%f ", m1[i]-m2[i]);
            if(m1[i]!=m2[i])
            {
                equals=false;
            }
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

    // Call Jacobi
    float* jacobiResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    jacobi(startVector, h, precomputedF, jacobiResult);
    end = get_wall_time();
    printf("Execution time Jacobi: %.3f seconds\n", end - start);

    // Call Gauss-Seidel
    float* gaussSeidelResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    gaussSeidel(startVector, h, precomputedF, gaussSeidelResult);
    end = get_wall_time();
    printf("Execution time Gauss-Seidel: %.3f seconds\n", end - start);

    // Call Gauss-Seidel Rot-Schwarz
    float* gaussSeidelRotSchwarzResult = malloc(size * size * sizeof(float));
    start = get_wall_time();
    gaussSeidelRotSchwarz(startVector, h, precomputedF, gaussSeidelRotSchwarzResult);
    end = get_wall_time();
    printf("Execution time Gauss-Seidel Rot-Schwarz: %.3f seconds", end - start);
    bool correct=true;
     correct=compare(gaussSeidelRotSchwarzResult,gaussSeidelResult);
        printf("is it correct: %s  \n" ,(correct)?"true":"false");

//Gaus Seidel Wavefront
    float* gaussSeidelWavefrontResult= malloc(size * size * sizeof(float));
    start = get_wall_time();
   gaussSeidelWavefront(startVector, h, precomputedF, gaussSeidelWavefrontResult);
    end = get_wall_time();
    printf("Execution time Gauss-Seidel Wavefront: %.3f seconds ", end - start);
    correct=compare(gaussSeidelWavefrontResult,gaussSeidelResult);
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
    printf("\nErgebnis analytisch:\n");
    printAnalyticalResult(h);*/

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
    for (i = 0; i < size * size; ++i)
    {
        a0[i] = startVector[i];
        a1[i] = startVector[i];
    }

    int k = 0;

//todo abbruchbedingung Max itera wieder ienführen
    for (k = 0; k < 1; k++)
    {
printf("2 \n");
        a0 = a1;
        int currentEle = 0;
        int border = 0;
        int durchlauf;
        for (durchlauf = 0; durchlauf<size+size-1 ; durchlauf++) //border > (size /2)
        {

            if (durchlauf > (size -1))
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
            printf("durchlauf %i ", durchlauf);
             printf("border %i ", border);
              printf("i%i ", i);
              printf("current ele %i ", currentEle);
            printf(" berechneter index%i \n",((durchlauf - border - i)  +( i+border )* size));

                a1[(durchlauf - border - i)  +( i + border)* size] = 0.25
                        * (a1[(durchlauf - border - i)  + (i + border - 1)* size]
                           + a1[(durchlauf - border - 1 - i)  + (i + border)* size]
                           + a0[(durchlauf - border + 1 - i)  + (i + border)* size]
                           + a0[(durchlauf - border - i)  + (i + border + 1)* size]
                           + h * h * functionTable[(durchlauf - border - i)  + (i + border)* size]);
            }
        }


    }
      printf("dödel \n");
        for (i = 0; i < size * size; ++i)
        {
            gaussSeidelResult[i] = a1[i];
        }
        free(a0);
        free(a1);

}
