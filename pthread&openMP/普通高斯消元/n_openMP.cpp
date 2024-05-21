#include <iostream>

#include <cstring>
#include <cstdlib>

#include <chrono>
#include <cstring>

#include<typeinfo>

#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>

using namespace std;
using namespace std::chrono;
const int n=1000;
const int pc=4;
float a[n][n];
void init(){
	memset(a,0,sizeof(a));
	for(int i=0;i<n;i++){
		a[i][i]=1;
		for(int j=i+1;j<n;j++)
			a[i][j]=rand();
	}
	for(int k=0;k<n;k++)
		for(int i=k+1;i<n;i++)
			for(int j=0;j<n;j++)
				a[i][j]+=a[k][j];
}

void gauss(){
	int i=0,j=0,k=0;
    float tmp=0;
    #pragma omp parallel  num_threads(pc),private(i,j,k,tmp)
    for (k=0;k<n;k++) {
        #pragma omp single
        {
			tmp=a[k][k];
        	for(j=k+1;j<n;j++){
            	a[k][j]=a[k][j]/tmp;
        	}
        	a[k][k]=1.0;
        }

        #pragma omp for schedule(static)
        for(i=k+1;i<n;i++){
            for(j=k+1;j<n;j++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
            a[i][k]=0;
        }
    }
}
int main(){
	init();
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	gauss();
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()<<std::endl;
	return 0;
}
