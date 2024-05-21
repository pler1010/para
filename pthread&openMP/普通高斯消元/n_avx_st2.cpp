#include <iostream>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <cstring>
#include<semaphore.h>
#include<pthread.h>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
using namespace std::chrono;
const int n=1000;
const int pc=8;
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
struct Param{
	int k,t_id;
};
pthread_barrier_t barrier_Division;
pthread_barrier_t barrier_Elimination;
void* threadf(void* param) {
    Param* p=(Param*)param;
    int t_id=p->t_id;

    for(int k=0;k<n;k++){ //0号线程做除法
        if(t_id==0){
            for(int j=k+1;j<n;j++){
                a[k][j]=a[k][j]/a[k][k];
            }
            a[k][k]=1.0;
        }

        pthread_barrier_wait(&barrier_Division);//第一个同步点

        for (int i=k+1+t_id;i<n;i+=pc) {
        	
            __m256 vaik = _mm256_set1_ps(a[i][k]);
            int j;
            for (j=k+1;j+8<=n;j+=8){
                __m256 vakj=_mm256_loadu_ps(&a[k][j]);
                __m256 vaij=_mm256_loadu_ps(&a[i][j]);
                __m256 vx=_mm256_mul_ps(vakj,vaik);
                vaij=_mm256_sub_ps(vaij,vx);
                _mm256_storeu_ps(&a[i][j],vaij);
            }
            for(;j<n;j++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
            a[i][k]=0;
        }

        pthread_barrier_wait(&barrier_Elimination);//第二个同步点

    }
    pthread_exit(NULL);
    return NULL;
}

void gauss() {
    pthread_barrier_init(&barrier_Division,NULL,pc);
    pthread_barrier_init(&barrier_Elimination,NULL,pc);

    pthread_t* handle = (pthread_t*)malloc(pc*sizeof(pthread_t));
    Param* param = (Param*)malloc(pc*sizeof(Param));

    for (int t_id=0;t_id<pc;t_id++) {
        param[t_id].t_id=t_id;
        param[t_id].k=0;
        pthread_create(&handle[t_id],NULL,threadf,&param[t_id]);

    }
    for (int t_id=0;t_id<pc;t_id++){
        pthread_join(handle[t_id],NULL);
    }

    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Elimination);

    free(handle);
    free(param);
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
