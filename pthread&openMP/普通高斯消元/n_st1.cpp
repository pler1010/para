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
sem_t sem_main;  //信号量
sem_t sem_workstart[pc];
sem_t sem_workend[pc];
struct Param{
	int k;
	int t_id;
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
void* threadFunc(void* param) {  //平凡，动态8线程+信号量
    Param* p=(Param*)param;
    int t_id=p->t_id;

    for (int k=0;k<n;k++){
        sem_wait(&sem_workstart[t_id]);//阻塞，等待主线程除法完成

        for (int i=k+1+t_id;i<n;i+=pc){
            for (int j=k+1;j<n;j++){
                a[i][j]=a[i][j]-a[i][k]*a[k][j];
            }
            a[i][k]=0;
        }

        sem_post(&sem_main);        //唤醒主线程
        sem_wait(&sem_workend[t_id]);  //阻塞，等待主线程唤醒进入下一轮

    }
    pthread_exit(NULL);
    return NULL;
}

void gauss(){
    sem_init(&sem_main,0,0); //初始化信号量
    for(int i=0;i<pc;i++) {
        sem_init(&sem_workend[i],0,0);
        sem_init(&sem_workstart[i],0,0);
    }
    pthread_t* handle=(pthread_t*)malloc(pc*sizeof(pthread_t));
    Param* param=(Param*)malloc(pc*sizeof(Param));
    for(int t_id=0;t_id<pc;t_id++) {
        param[t_id].t_id=t_id;
        param[t_id].k=0;
        pthread_create(&handle[t_id],NULL,threadFunc,&param[t_id]);
    }

    for(int k=0;k<n;k++){
        for(int j=k+1;j<n;j++){
            a[k][j]=a[k][j]/a[k][k];
        }
        a[k][k]=1.0;

        for(int t_id=0;t_id<pc;t_id++){  //唤起子线程
            sem_post(&sem_workstart[t_id]);
        }

        for (int t_id=0;t_id<pc;t_id++){  //主线程睡眠
            sem_wait(&sem_main);
        }

        for (int t_id=0;t_id<pc;t_id++){  //再次唤起工作线程，进入下一轮消去
            sem_post(&sem_workend[t_id]);
        }

    }
    for (int t_id=0;t_id<pc;t_id++){
        pthread_join(handle[t_id],NULL);
    }
    sem_destroy(&sem_main);    //销毁线程
    for (int t_id=0;t_id<pc;t_id++)
        sem_destroy(&sem_workstart[t_id]);
    for (int t_id=0;t_id<pc;t_id++)
        sem_destroy(&sem_workend[t_id]);

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
