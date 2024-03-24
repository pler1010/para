#include <iostream>
#include <ratio>
#include <chrono>
#include <cstring>
using namespace std::chrono;
#define uint unsigned int
#define add_1(x) sum+=a[x];sum+=a[x|1];
#define add_2(x) add_1(x);add_1(x|2);
#define add_3(x) add_2(x);add_2(x|4);
#define add_4(x) add_3(x);add_3(x|8);
#define add_5(x) add_4(x);add_4(x|16);
#define add_6(x) add_5(x);add_5(x|32);
#define add_7(x) add_6(x);add_6(x|64);
#define add_8(x) add_7(x);add_7(x|128);
#define add_9(x) add_8(x);add_8(x|256);
#define add_10(x) add_9(x);add_9(x|512);
#define add_11(x) add_10(x);add_10(x|1024);
#define add_12(x) add_11(x);add_11(x|2048);
#define add_13(x) add_12(x);add_12(x|4096);
#define add_14(x) add_13(x);add_13(x|8192);
#define add_15(x) add_14(x);add_14(x|16384);
#define add_16(x) add_15(x);add_15(x|32768);
#define add_17(x) add_16(x);add_16(x|65536);
#define add_18(x) add_17(x);add_17(x|131072);
uint add6(int n,uint a[]){
	uint sum=0;
	add_12(0);
	return sum;
}
const int maxn=1e7+10;
uint a[maxn];
uint b[maxn];
void init(int n,uint a[]){
	memcpy(a,b,n*4);
}
int main(){
	int n=1<<12;
	int T=1000000000ll/n;
	for(int i=0;i<n;i++) b[i]=i*i+i+10;
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	for(int i=1;i<=T;i++) init(n,a),add6(n,a);
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()/T/n*1000000000<<std::endl;
	return 0;
}
