#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <vector>
#include<pthread.h>
using namespace std;
using namespace std::chrono;
#define N 130
const int ndn=17;
const int maxn=17;
#define mm 22
#define nn 8

#define path1 "./../../Groebner/1/nd.txt"

#define path2 "./../../Groebner/1/wait.txt"

int geth(){
	return mm;
}
int getx(){
	return nn;
}
std::vector<int> b[N+nn];
void getin(vector<int> &x){
	int tmp=0;
	char ch=getchar();
	while(ch<'0'||ch>'9') ch=getchar();
	while((ch>='0'&&ch<='9')||ch==' '){
		if(ch==' '){
			x.push_back(tmp);
			tmp=0;
		}
		else{
			tmp=tmp*10+ch-'0';
		}
		ch=getchar();
	}
}
void print(const vector<int> &x){
	for(int i=0;i<x.size();i++) printf("%d ",x[i]);
	puts("");
}
const int pc=8;
struct Param{
	int rank;
	int st;
	int to;
	vector<int> *ve1;
	const vector<int> *ve2;
};
void* threadFunc(void *param){
	Param *p=(Param *)param;
	vector<int> tmp;
	vector<int> &x=*(p->ve1);
	const vector<int> &y=*(p->ve2);
	auto be1=x.begin()+p->st;
	auto en1=x.begin()+p->to;
	auto be2=y.end();
	if(be1!=x.end()) be2=lower_bound(y.begin(),y.end(),*be1);
	auto en2=y.end();
	if(en1!=x.end()) en2=lower_bound(y.begin(),y.end(),*en1);
	auto it1=be1;
	auto it2=be2;
	while(it1<en1||it2<en2){
		if(it1==en1||it2==en2){
			if(it1==en1){
				tmp.push_back(*it2);
				it2++;
			}
			else{
				tmp.push_back(*it1);
				it1++;
			}
			continue;
		}
		if(*it1==*it2){
			it1++;
			it2++;
			continue;
		}
		else if(*it1<*it2){
			tmp.push_back(*it2);
			it2++;
		}
		else{
			tmp.push_back(*it1);
			it1++;
		}
	}
	for(int i=0;i<tmp.size();i++) x.push_back(tmp[i]);
	return NULL;
}
void insert(vector<int> &x){
	for(int i=N-1;i>=0;i--){
		if(x[0]==i){
			if(b[i].size()){
				pthread_t* handle=(pthread_t*)malloc(pc*sizeof(pthread_t));
				Param* param=(Param*)malloc(pc*sizeof(Param));
				for(int id=0;id<pc;id++){
					param[id].st=N/pc*id;
					param[id].to=N/pc*(id+1);
					if(id==pc-1) param[id].to=N;
					param[id].ve1=&x;
					param[id].ve2=&b[i];
				};
				for(int id=0;id<pc;id++){
					pthread_create(&handle[id],NULL,threadFunc,&param[id]);
				}
				for(int id=0;id<pc;id++) pthread_join(handle[id],NULL);
				free(handle);
				free(param);
			}
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
