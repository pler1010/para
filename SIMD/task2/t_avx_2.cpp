#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>

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
#define N 130
const int ndn=17;
const int maxn=24;
#define mm 22
#define nn 8

#define path1 "Groebner/1-130-22-8/nd.txt"

#define path2 "Groebner/1-130-22-8/wait.txt"

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
	void calcxor(Bitset &b){
		int i=0;
		while(i<ndn&&(((int)(&a[i]))&31)){
			a[i]^=b.a[i];
			i++;
		}
		for(;i+7<ndn;i+=8){ 
			__m256i t1=_mm256_load_si256((__m256i*)(&a[i]));
			__m256i t2=_mm256_load_si256((__m256i*)(&b.a[i]));
			t1=_mm256_xor_si256(t1,t2);
			_mm256_store_si256((__m256i*)(&a[i]),t1);
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
void insert(Bitset &x){
	for(int i=N-1;i>=0;i--){
		if(x.get(i)){
			if(b[i].get(i)) x.calcxor(b[i]);
			else{
				b[i]=x;
				break;
			}
		}
	}
}
void solve(){
	int m=geth(),n=getx();
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
	for(int i=0;i<n;i++) insert(b[N+i]);
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
