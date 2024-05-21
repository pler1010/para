#include <iostream>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <cstring>
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

void* threadFunc(void* param) {
    Param* p = (Param*)param;
    int k=p->k;           //消去的轮次
    int t_id=p->t_id;     //线程
    int i=k+t_id+1;   //获取任务
    for (;i<n;i+=pc){
        for(int j=k+1;j<n;j++){
            a[i][j]=a[i][j]-a[i][k]*a[k][j];
        }
        a[i][k]=0;
    }
    pthread_exit(NULL);
    return NULL;
}

void gauss() {    //串行，动态分配8线程
    for(int k=0;k<n;k++){
        for(int j=k+1;j<n;j++) a[k][j]=a[k][j]/a[k][k];
        a[k][k]=1.0;
        
        pthread_t* handle=(pthread_t*)malloc(pc*sizeof(pthread_t));
        Param* param=(Param*)malloc(pc*sizeof(Param));

        for(int t_id=0;t_id<pc;t_id++){//分配任务
            param[t_id].k=k;
            param[t_id].t_id=t_id;
        }

        for (int t_id=0;t_id<pc;t_id++){
            pthread_create(&handle[t_id],NULL,threadFunc,&param[t_id]);
        }

        for (int t_id=0;t_id<pc;t_id++) {
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
