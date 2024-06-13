#include <iostream>

#include <cstring>
#include <cstdlib>

#include <chrono>
#include <cstring>

#include<typeinfo>

#include <sycl/sycl.hpp>

using namespace std;
using namespace std::chrono;
using namespace sycl;
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
void gauss_buffer(){
	std::vector<float> aa(n*n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			aa[i*n+j]=a[i][j];
	buffer<float,2> buf(aa.data(),range<2>(n,n));
	queue q{property::queue::in_order()};
	device my_device = q.get_device();
	std::cout << "Device: " << my_device.get_info<info::device::name>() << std::endl;
 
	int n = buf.get_range()[0];
	for (int k = 0; k < n; k++) {
 
		q.submit([&](handler& h) {
			accessor m{ buf, h, read_write };
			h.parallel_for(range(n - k), [=](auto idx) {
				int j = k + idx;
				m[k][j] = m[k][j] / m[k][k];
			});
		});
 
		q.submit([&](handler& h) {
			accessor m{ buf, h, read_write };
			h.parallel_for(range(n - (k + 1), n - (k + 1)), [=](auto idx) {
				int i = k + 1 + idx.get_id(0);
				int j = k + 1 + idx.get_id(1);
				m[i][j] = m[i][j] - m[i][k] * m[k][j];
			});
		});
 
		q.submit([&](handler& h) {
			accessor m{ buf, h, read_write };
			h.parallel_for(range(n - (k + 1)), [=](auto idx) {
				int i = k + 1 + idx;
				m[i][k] = 0;
			});
		});
	}
	q.wait();
	host_accessor h_a(buf,read_only);
	for(int i=0;i<n;i++) for(int j=0;j<n;j++) a[i][j]=h_a[i][j];
}
void gauss_usm(){
	queue q{property::queue::in_order()};
	device my_device = q.get_device();
	std::cout << "Device: " << my_device.get_info<info::device::name>() << std::endl;
 
 	float *m=malloc_device<float>(n*n,q);
 	q.memcpy(m,a,sizeof(float)*n*n).wait();
	for (int k = 0; k < n; k++) {
 
		q.submit([&](handler& h) {
			h.parallel_for(range(n - k), [=](auto idx) {
				int j = k + idx;
				m[k*n+j] = m[k * n + j] / m[k * n + k];
			});
		});
 
		q.submit([&](handler& h) {
			h.parallel_for(range(n - (k + 1), n - (k + 1)), [=](auto idx) {
				int i = k + 1 + idx.get_id(0);
				int j = k + 1 + idx.get_id(1);
				m[i * n + j] = m[i * n + j] - m[i * n + k] * m[k * n + j];
			});
		});
 
		q.submit([&](handler& h) {
			h.parallel_for(range(n - (k + 1)), [=](auto idx) {
				int i = k + 1 + idx;
				m[i * n + k] = 0;
			});
		});
	}
	q.wait();
 	q.memcpy(a,m,sizeof(float)*n*n).wait();
}
int main(){
	init();
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	gauss();
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()<<std::endl;
	
    init();
    gauss_buffer();
	
    init();
	t1=high_resolution_clock::now();
	gauss_buffer();
	t2=high_resolution_clock::now();
	time_span=t2-t1;
	std::cout<<time_span.count()<<std::endl;
	
	init();
	t1=high_resolution_clock::now();
	gauss_usm();
	t2=high_resolution_clock::now();
	time_span=t2-t1;
	std::cout<<time_span.count()<<std::endl;
	return 0;
}

