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
const int pc=8;
float a[n][n];

struct Param{
	int k,t_id;
};

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

void* threadFunc(void* param) {  //动态分配8线程
    Param* p=(Param*)param;
    int k=p->k;           //消去的轮次
    int t_id=p->t_id;     //线程
    int i=k+t_id+1;   //获取任务

    for (;i<n;i+=pc) {
        __m256 vaik=_mm256_set1_ps(a[i][k]);
        int j;
        for (j=k+1;j+8<=n;j+=8){
            __m256 vakj=_mm256_loadu_ps(&a[k][j]);
            __m256 vaij=_mm256_loadu_ps(&a[i][j]);
            __m256 vx=_mm256_mul_ps(vakj,vaik);
            vaij=_mm256_sub_ps(vaij,vx);
            _mm256_storeu_ps(&a[i][j],vaij);
        }
        for (;j<n;j++) a[i][j]=a[i][j]-a[i][k]*a[k][j];
        a[i][k] = 0;
    }
    pthread_exit(NULL);
    return NULL;
}

void gauss() {            //neon,动态分配8线程
    for (int k=0;k<n;k++) {
        __m256 vt=_mm256_set1_ps(a[k][k]);
        int j=0;
        for (j=k+1;j+8<=n;j+=8){
            __m256 va=_mm256_loadu_ps(&a[k][j]);
            va=_mm256_div_ps(va,vt);
            _mm256_storeu_ps(&a[k][j],va);
        }
        for (;j<n;j++) a[k][j]=a[k][j]/a[k][k];
        a[k][k] = 1.0;

        pthread_t* handle=(pthread_t*)malloc(pc*sizeof(pthread_t));
        Param* param=(Param*)malloc(pc*sizeof(Param));

        for(int t_id=0;t_id<pc;t_id++){//分配任务
            param[t_id].k=k;
            param[t_id].t_id=t_id;
        }

        for(int t_id=0;t_id<pc;t_id++){
            pthread_create(&handle[t_id],NULL,threadFunc,&param[t_id]);
        }

        for(int t_id=0;t_id<pc;t_id++){
            pthread_join(handle[t_id],NULL);
        }
        free(handle);
        free(param);
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
