#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <vector>
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
/*class Bitset{
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
		for(int i=0;i<ndn;i++) a[i]^=b.a[i];
	}
};
Bitset b[N+nn];*/
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
void calcxor(vector<int> &x,const vector<int> &y){
	vector<int> tmp;
	auto it1=x.begin();
	auto it2=y.begin();
	while(it1!=x.end()||it2!=y.end()){
		if(it1==x.end()||it2==y.end()){
			if(it1==x.end()){
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
	x=tmp;
}
void insert(vector<int> &x){
	for(int i=N-1;i>=0;i--){
		if(x[0]==i){
			if(b[i].size()) calcxor(x,b[i]);
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
