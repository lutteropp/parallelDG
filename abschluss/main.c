#include <stdio.h>
#include <stdlib.h>

float f1(float x, float y);
float* jacobi(float* startVector, float h, float* f);
float* gaussSeidel(float * startVector, float h, float* f);
float* computeFunctionTable(float h);
void printResultMatrix(float* matrix);
void printAnalyticalResult(float h);

int size;

#define CO(i,j) ( (j) * (size) + (i) )

float* jacobi(float* startVector, float h, float* f) {
	float* a0 = startVector; // last iteration
	float* a1 = startVector; // current iteration
	int i, j, k;
	//todo abbruchbedingung
	for (k = 0; k < 1000000; ++k) {
		a0 = a1;
		for (j = 1; j < size; j++) {
			for (i = 1; i < size; i++) {
				a1[CO(i,j)] = a0[CO(i, j - 1)]
						+ a0[CO(i - 1, j)]
						+ a0[CO(i, j + 1)]
						+ a0[CO(i + 1, j)]
						+ h*h * f[CO(i, j)];
				a1[CO(i,j)] *= 0.25;
			}
		}
	}
	return a1;
}

float* gaussSeidel(float * startVector, float h, float* f) {
	float* a0 = startVector; // last iteration
	float* a1 = startVector; // current iteration
	int i, j, k;
	//todo abbruchbedingung
	for (k = 0; k < 1000000; ++k) {
		a0 = a1;
		for (j = 1; j < size; j++) {
			for (i = 1; i < size; i++) {
				a1[CO(i,j)] = a1[CO(i, j - 1)]
						+ a1[CO(i - 1, j)]
						+ a0[CO(i, j + 1)]
						+ a0[CO(i + 1, j)]
						+ h*h * f[CO(i, j)];
				a1[CO(i,j)] *= 0.25;
			}
		}
	}
	return a1;
}


float* computeFunctionTable(float h) {
	float* functionTable = (float*) malloc(size * size* sizeof(float));
	int i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			float x1 = i * h;
			float y1 = j * h;
			functionTable[CO(i, j)] = f1(x1, y1);
		}
	}
	return functionTable;
}

float f1(float x, float y) {
	return (32 * (x * (1 - x) + y * (1 - y)));
}

void printResultMatrix(float* matrix) {
	int i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			float val = matrix[CO(i, j)];
			if (j == 0) {
				printf("%.3f", val);
			} else if (j < size - 1) {
				printf(" %.3f ", val);
			} else {
				printf(" %.3f\n", val);
			}
		}
	}
}

void printAnalyticalResult(float h) {
	int i, j;
	for (i = 0; i < size; ++i) {
		float x = i * h;
		for (j = 0; j < size; ++j) {
			float y = j * h;
			float val;
			if (i == 0 || j == 0 || i == size - 1 || j == size - 1) {
				val = 0; // Randbedingung
			} else {
				val = 16 * x * (1 - x) * y * (1 - y);
			}

			if (j == 0) {
				printf("%.3f", val);
			} else if (j < size - 1) {
				printf(" %.3f ", val);
			} else {
				printf(" %.3f\n", val);
			}
		}
	}
}

int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Usage: ./parallelDG SIZE\n (SIZE stands for the number of grid points in each dimension)\n");
		return -1;
	}

	size = atoi(argv[1]);
	if (size <= 1) {
		printf("Error! Size must be greater than 1.\n");
		return -1;
	}
	float h = 1 / (float) (size - 1);
	// Precompute function f
	float* precomputedF = computeFunctionTable(h);
	// Create random start vector with zeroes at the proper positions (at the border)
	float* startVector = (float*) malloc(size * size * sizeof(float));
	srand(time(NULL));
	int i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			float val = 0;
			if (i != 0 && j != 0) {
				val = (float) rand() / RAND_MAX;
			}
			startVector[CO(i, j)] = val;
		}
	}
	// Call Jacobi
	float* jacobiResult = jacobi(startVector, h, precomputedF);
	// Call Gauss-Seidel
	float* gaussSeidelResult = gaussSeidel(startVector, h, precomputedF);

	// TODO: Check for correctness and compare both results

	// TODO: The following is just debug code. Remove afterwards.
	printf("Ergebnis Jacobi-Verfahren:\n");
	printResultMatrix(jacobiResult);
	printf("\nErgebnis Gauss-Seidel-Verfahren:\n");
	printResultMatrix(gaussSeidelResult);
	printf("\nErgebnis analytisch:\n");
	printAnalyticalResult(h);


	free(jacobiResult);
	free(gaussSeidelResult);
	free(startVector);
	free(precomputedF);
	return 0;
}

float* wavefront(float * matrix, float h, float* f) {
  float* a0 = matrix;
  float* a1 = matrix;
  int k = 0;
  int size = 1 / h;
//todo abbruchbedingung
  for (k = 0; k < 100000; k++)
  {

    a0 = a1;
    int currentEle = 1;
    int border = 0;
    int durchlauf;
    for (durchlauf = 0; border > (size / 2); durchlauf++)
    {
      if (currentEle >= (size / 2))
      {
        currentEle--;
        border++;
      } else
      {
        currentEle++;
      }
      int i = 0;
      for (i = 0; i < currentEle; i++)
      {
        a1[(durchlauf + border - i) * size + i + border] = 0.25
            * (a1[(durchlauf + border - i) * size + i + border - 1]
                + a1[(durchlauf + border - 1 - i) * size + i + border]
                + a0[(durchlauf + border + 1 - i) * size + i + border]
                + a0[(durchlauf + border - i) * size + i + border + 1]
                + h * h * f[(durchlauf + border - i) * size + i + border]);
      } //indexe machen
    }

//todo erinnern was x war
    // a1[i*size + j]=0.25*(a1[i*size + (j-1)] + a1[(i-1)*size + j]+a0[(i+1)*size + j]+a0[i*size + (j+1)]+h*h*f[i*size+j]);

  }
  return a1;
}
