#include <iostream>
#include <cstring>
#include <chrono>
#include <ratio>
using namespace std::chrono;
const int maxn=1e4+10;
#define uint unsigned int 
uint sum[maxn];
uint b[maxn][maxn];
uint a[maxn];
uint tb[maxn][maxn];
uint ta[maxn];
void init(int n){
	memcpy(a,ta,n*4);
	for(int i=0;i<n;i++)
		memcpy(b[i],tb[i],n*4);
}
void mul0(int n){
	for(int i=0;i<n;i++){
		sum[i]=0;
		for(int j=0;j<n;j++) sum[i]+=b[j][i]*a[j];
	}
}
void mul1(int n){
	for(int i=0;i<n;i++) sum[i]=0;
	for(int j=0;j<n;j++)
		for(int i=0;i<n;i++)
			sum[i]+=b[j][i]*a[j];
}
#define solve1_0(x) sum[i]+=b[x][i]*a[x];
#define solve1_1(x) solve1_0(x);solve1_0(x|1);
#define solve1_2(x) solve1_1(x);solve1_1(x|2);
#define solve1_3(x) solve1_2(x);solve1_2(x|4);
#define solve1_4(x) solve1_3(x);solve1_3(x|8);
#define solve1_5(x) solve1_4(x);solve1_4(x|16);
#define solve1_6(x) solve1_5(x);solve1_5(x|32);
#define solve1_7(x) solve1_6(x);solve1_6(x|64);
#define solve1_8(x) solve1_7(x);solve1_7(x|128);
#define solve1_9(x) solve1_8(x);solve1_8(x|256);
#define solve1_10(x) solve1_9(x);solve1_9(x|512);
#define solve1_11(x) solve1_10(x);solve1_10(x|1024);
#define solve1_12(x) solve1_11(x);solve1_11(x|2048);

#define solve2_0(x) sum[x]+=b[j][x]*a[j];
#define solve2_1(x) solve2_0(x);solve2_0(x|1);
#define solve2_2(x) solve2_1(x);solve2_1(x|2);
#define solve2_3(x) solve2_2(x);solve2_2(x|4);
#define solve2_4(x) solve2_3(x);solve2_3(x|8);
#define solve2_5(x) solve2_4(x);solve2_4(x|16);
#define solve2_6(x) solve2_5(x);solve2_5(x|32);
#define solve2_7(x) solve2_6(x);solve2_6(x|64);
#define solve2_8(x) solve2_7(x);solve2_7(x|128);
#define solve2_9(x) solve2_8(x);solve2_8(x|256);
#define solve2_10(x) solve2_9(x);solve2_9(x|512);
#define solve2_11(x) solve2_10(x);solve2_10(x|1024);
#define solve2_12(x) solve2_11(x);solve2_11(x|2048);
void mul2(int n){
	for(int i=0;i<n;i++){
		sum[i]=0;
		solve1_12(0);
		//for(int j=0;j<n;j++) sum[i]+=b[j][i]*a[j];
	}
}
void mul3(int n){
	for(int i=0;i<n;i++) sum[i]=0;
	for(int j=0;j<n;j++){
		solve2_12(0);
		//for(int i=0;i<n;i++){
		//	sum[i]+=b[j][i]*a[j];
		//}
	}
}
int main(){
	int n=1<<12;
	int T=1000000000ll/(n*n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++)
			tb[i][j]=i+j;
		ta[i]=i;
	}
	init(n);
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	for(int i=1;i<=T;i++) mul2(n);
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()/T/n/n*1000000000<<std::endl;
	return 0;
}
