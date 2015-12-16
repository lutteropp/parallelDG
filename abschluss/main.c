#include <stdio.h>
#include <stdlib.h>
float f1(float x,float y);
float* createTestMatrix(int h);
float* jacobi(float * matrix,float* u0,float h, float* x,float*y);
float f1(float x,float y);

int main()
{
    printf("Hello world!asdf\n");
    return 0;
}
float* createTestMatrix(int h)
{
    int size = 1/h;
    float* a = (float *)malloc(size *size * sizeof(float));
    int i;
    for(i =0;  i<size; i++)
    {
        int j;
        for(j =0; j<size; j++)
        {
            if(i==j)
                a[i*size + j] =4;
        }
        if (i==(j-1)||(i-1)==j)
        {
            a[i*size + j] =1;
        }

    }
    return a;

}
float* jacobi(float * matrix,float* u0,float h, float* x,float*y)
{
    float* a0=matrix;
    float* a1=matrix;
    int k =0;
    int size=1/h;
    //todo abbruchbedingung
    for( k =0; k<100000; k++)
    {
        a0=a1;
        int j;
        for ( j=1; j<size; j++)
        {
            int i;
            for( i=1; i<size; i++)
            {
                //todo erinnern was x war
                a1[i*size + j]=0,25*(a0[i*size + (j-1)] + a0[(i-1)*size + j]+a0[(i+1)*size + j]+a0[i*size + (j+1)]+h*h*f1(x[i*size + j],y[i*size + j]));
            }
        }
    }
    return a1;
}
float f1(float x,float y)
{
    float asdf=(32*(x*(1-x)+y*(1-y)));
    return asdf;
}
