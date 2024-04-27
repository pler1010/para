#include <iostream>

#include <cstring>
#include <cstdlib>

#include <chrono>
#include <cstring>

#include<arm_neon.h>

using namespace std;
using namespace std::chrono;
const int n=1000;
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
	for(int i=0;i<n;i++){
		float32x4_t t1=vdupq_n_f32(a[i][i]);
		int j=0;
		for(j=i+1;j+3<n;j+=4){
			float32x4_t t2=vld1q_f32(&a[i][j]);
			t2=vdivq_f32(t2,t1);
			vst1q_f32(&a[i][j],t2);
		}
		for(;j<n;j++) a[i][j]=a[i][j]/a[i][i];
		a[i][i]=1;
		for(int j=i+1;j<n;j++){
			float32x4_t tmp=vdupq_n_f32(a[j][i]);
			int k=0;
			for(k=i+1;k+3<n;k+=4){
				float32x4_t tmpi=vld1q_f32(&a[i][k]);
				float32x4_t tmpj=vld1q_f32(&a[j][k]);
				float32x4_t ntmp=vmulq_f32(tmp,tmpi);
				tmpj=vsubq_f32(tmpj,ntmp);
				vst1q_f32(&a[j][k],tmpj);
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
