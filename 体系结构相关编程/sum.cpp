#include <iostream>
#include <ratio>
#include <chrono>
#include <cstring>
using namespace std::chrono;
#define uint unsigned int
uint add0(int n,uint a[]){
	uint sum=0;
	for(int i=0;i<n;i++) sum+=a[i];
	return sum;
}
uint add1(int n,uint a[]){
	uint sum1=0,sum2=0;
	for(int i=0;i<n;i+=2){
		sum1+=a[i];
		sum2+=a[i+1];
	}
	return sum1+sum2;
}
void recursion(int n,uint a[]){
	if(n==1) return;
	for(int i=0;i<n/2;i++) a[i]+=a[n-i-1];
	recursion(n/2,a);
}
uint add2(int n,uint a[]){
	recursion(n,a);
	return a[0];
}
uint add3(int n,uint a[]){
	for(int m=n;m>1;m/=2){
		for(int i=0;i<m/2;i++){
			a[i]=a[i*2]+a[i*2+1];
		}
	}
	return a[0];
}
uint add4(int n,uint a[]){
	uint sum=0;
	for(int i=0;i<n;i+=2) sum+=a[i],sum+=a[i+1];
	return sum;
}
uint add5(int n,uint a[]){
	for(int m=n;m>1;m/=2){
		for(int i=0;i<m/2;i++){
			a[i]+=a[i+m/2];
		}
	}
	return a[0];
}
const int maxn=1e7+10;
uint a[maxn];
uint b[maxn];
void init(int n,uint a[]){
	memcpy(a,b,n*4);
}
int main(){
	int n=1<<18;
	int T=1000000000ll/n;
	for(int i=0;i<n;i++) b[i]=i*i+i+10;
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	for(int i=1;i<=T;i++) init(n,a),add0(n,a);
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()/T/n*1000000000<<std::endl;
	return 0;
}
