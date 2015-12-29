#include <stdio.h>
#include <stdlib.h>
float f1(float x,float y);
float* createTestMatrix(int h);
float* jacobi(float * matrix,float h, float* f);
float* gausseidel(float * matrix,float h, float* f);
float* f(float*x,float*y,float h);


int main()
{
printf("Hello world!asdf\n");
    float h = 2;
    float* otto= createTestMatrix(h);

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
float* gausseidel(float * matrix,float h, float* f){

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
                a1[i*size + j]=0,25*(a1[i*size + (j-1)] + a1[(i-1)*size + j]+a0[(i+1)*size + j]+a0[i*size + (j+1)]+h*h*f[i*size+j]);
            }
        }
    }
    return a1;
}

float* jacobi(float * matrix,float h, float* f)
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
                a1[i*size + j]=0,25*(a0[i*size + (j-1)] + a0[(i-1)*size + j]+a0[(i+1)*size + j]+a0[i*size + (j+1)]+h*h*f[i*size+j]);
            }
        }
    }
    return a1;
}
float* f(float*x,float*y,float h)
{
    int size= 1/h;
    float* tmp=x;
    int i=0;

    for( i=0; i<size; i++)
    {
        int j=0;
        for(j =0; j<size; j++)
        {
        float x1=x[i*size+j];
        float y1= y[i*size+j];
            tmp[i*size+j]=f1(x1,y1);
        }
    }
    return tmp;
}
float f1(float x,float y)
{
    float asdf=(32*(x*(1-x)+y*(1-y)));
    return asdf;
}
/*
float* wavefront(float * matrix,float h, float* f)
{
    float* a0=matrix;
    float* a1=matrix;
    int k =0;
    int size=1/h;
    //todo abbruchbedingung
    for( k =0; k<100000; k++)
    {
        
        a0=a1;
        int currentEle=1;
        int border =0;
        int durchlauf;
        for(durchlauf =0;border>(n/2);durchlauf++){
            if(currentEle>=(n/2)){
                currentEle--;
                border++;
            }else{
                currentEle++;
            }
            int i =0;
            for(i=0;i<currentEle;i++){
                a1[(durchlauf+border-i*size + i+border]=0,25*(a1[i*size + (j-1)] + a1[(i-1)*size + j]+a0[(i+1)*size + j]+a0[i*size + (j+1)]+h*h*f[i*size+j]);
            }//indexe machen
        }
        
                //todo erinnern was x war
                 a1[i*size + j]=0,25*(a1[i*size + (j-1)] + a1[(i-1)*size + j]+a0[(i+1)*size + j]+a0[i*size + (j+1)]+h*h*f[i*size+j]);
       
    }
    return a1;
}

current elemente =1
border =0;
for durchauf ++ bis border >n/2

if (current ele >=n/2) curre ele--; border++; else cur el ++;
 for(current elemente auf der diagonalen
	int lala = (durchlauf -i)(i)
 A(durchlauf+border -i)(i+border) =bla blub


*/
