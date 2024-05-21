#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
#include<pthread.h>
using namespace std;
using namespace std::chrono;
#define N 130
const int ndn=17;
const int maxn=17;
#define mm 22
#define nn 8

const int pc=8;

int n,m;

#define path1 "./../../Groebner/1/nd.txt"

#define path2 "./../../Groebner/1/wait.txt"

int geth(){
	return mm;
}
int getx(){
	return nn;
}
class Bitset{
	unsigned int a[maxn];
public:
	Bitset(){
		memset(a,0,sizeof(a));
	}
	void init(){
		memset(a,0,sizeof(a));
	}
	void set(int pst,bool value){
		if(get(pst)!=value) a[pst/32]^=1u<<(pst&31);
	}
	bool get(int pst) const {
		return (a[pst/32]>>(pst&31))&1;
	}
	void calcxor(const Bitset &b){
		int i;
		for(i=0;i+7<ndn;i+=8){ 
			__m256i t1=_mm256_loadu_si256((__m256i*)(&a[i]));
			__m256i t2=_mm256_loadu_si256((__m256i*)(&b.a[i]));
			t1=_mm256_xor_si256(t1,t2);
			_mm256_storeu_si256((__m256i*)(&a[i]),t1);
		}
		for(;i<ndn;i++) a[i]^=b.a[i];
	}
};
Bitset b[N+nn];
void getin(Bitset &x){
	x.init();
	int tmp=0;
	char ch=getchar();
	while(ch<'0'||ch>'9') ch=getchar();
	while((ch>='0'&&ch<='9')||ch==' '){
		if(ch==' '){
			x.set(tmp,1);
			tmp=0;
		}
		else{
			tmp=tmp*10+ch-'0';
		}
		ch=getchar();
	}
}
void print(const Bitset &x){
	for(int i=N-1;i>=0;i--)
		if(x.get(i))
			printf("%d ",i);
	puts("");
}
pthread_mutex_t latch[N];
void insert(Bitset &x){
	for(int i=N-1;i>=0;i--){
		pthread_mutex_lock(&latch[i]);
		if(x.get(i)){
			if(b[i].get(i)) x.calcxor(b[i]);
			else{
				b[i]=x;
				pthread_mutex_unlock(&latch[i]);
				break;
			}
		}
		pthread_mutex_unlock(&latch[i]);
	}
}
struct Param{
	int t_id;
};
void* thread(void *p){
	for(int i=((Param *)p)->t_id;i<n;i+=pc) insert(b[N+i]);
	return NULL;
}
void start(){
    pthread_t* handle=(pthread_t*)malloc(pc*sizeof(pthread_t));
    Param* param=(Param*)malloc(pc*sizeof(Param));
    for(int i=0;i<N;i++) pthread_mutex_init(&latch[i],NULL);
    for (int t_id=0;t_id<pc;t_id++) {
        param[t_id].t_id=t_id;
        pthread_create(&handle[t_id],NULL,thread,&param[t_id]);

    }
	for (int t_id=0;t_id<pc;t_id++) {
        pthread_join(handle[t_id], NULL);
    }
    for(int i=0;i<N;i++) pthread_mutex_destroy(&latch[i]);
    free(handle);
    free(param);
}
void solve(){
	m=geth(),n=getx();
	freopen(path1,"r",stdin);
	for(int i=0;i<m;i++){
		getin(b[N]);
		insert(b[N]);
	}
	fclose(stdin);
	freopen(path2,"r",stdin);
	freopen("r.txt","w",stdout);
	for(int i=0;i<n;i++) getin(b[N+i]);
	
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	start();
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	cerr<<time_span.count()<<endl;
	
	for(int i=0;i<n;i++) print(b[N+i]);
	fclose(stdin);
	fclose(stdout);
}
int main(){
	solve();
	return 0;
}
