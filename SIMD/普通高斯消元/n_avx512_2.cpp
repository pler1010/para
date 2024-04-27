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
const int sn=1024;
float a[sn][sn];
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
	for(int i=0;i<n;i++){
		__m512 t1=_mm512_set1_ps(a[i][i]);
		int j=i+1;
		while(j<n&&(int(&a[i][j])&63)){
			a[i][j]=a[i][j]/a[i][i];
			j++;
		}
		for(;j+15<n;j+=16){
			__m512 t2=_mm512_load_ps(&a[i][j]);
			t2=_mm512_div_ps(t2,t1);
			_mm512_store_ps(&a[i][j],t2);
		}
		for(;j<n;j++) a[i][j]=a[i][j]/a[i][i];
		a[i][i]=1;
		for(int j=i+1;j<n;j++){
			__m512 tmp=_mm512_set1_ps(a[j][i]);
			int k=i+1;
			while(k<n&&(int(&a[j][k])&63)){
				a[j][k]-=a[j][i]*a[i][k];
				k++;
			}			
			for(;k+15<n;k+=16){
				__m512 tmpi=_mm512_load_ps(&a[i][k]);
				__m512 tmpj=_mm512_load_ps(&a[j][k]);
				__m512 ntmp=_mm512_mul_ps(tmp,tmpi);
				tmpj=_mm512_sub_ps(tmpj,ntmp);
				_mm512_store_ps(&a[j][k],tmpj);
			}
			for(;k<n;k++) a[j][k]-=a[j][i]*a[i][k];
			a[j][i]=0;
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
