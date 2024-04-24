#include <iostream>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <cstring>
using namespace std::chrono;
const int n=1000;
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
void gauss(){
	for(int i=0;i<n;i++){
		float tmp=a[i][i];
		for(int j=i;j<n;j++){
			a[i][j]/=tmp;
		}
		for(int j=i+1;j<n;j++){
			float tmp=a[j][i];
			for(int k=i;k<n;k++){
				a[j][k]-=tmp*a[i][k];
			}
		}
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
