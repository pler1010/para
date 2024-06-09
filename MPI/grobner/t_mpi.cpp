#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<vector>
#include<windows.h>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
#include<mpi.h>
#include<omp.h>
using namespace std;

#define NUM_THREADS 8


const int maxsize = 3000;
const int maxrow = 40000;
const int numBasis = 40000;
int num;

vector<int> tmpAns;

long long head, tail, freq;
map<int, int*>ans;			//答案

fstream RowFile("wait.txt", ios::in | ios::out);
fstream BasisFile("nd.txt", ios::in | ios::out);

ofstream out_mpi("res.txt");

int gRows[maxrow][maxsize];   //被消元行最多60000行，3000列
int gBasis[numBasis][maxsize];  //消元子最多40000行，3000列
int answers[maxrow][maxsize]; //存储消元完毕的行
map<int, int>firstToRow; //记录answers的每行和首项的对应关系

int ifBasis[numBasis] = { 0 };
int ifDone[maxrow] = { 0 };

void reset() {
	//	read = 0;
	memset(gRows, 0, sizeof(gRows));
	memset(gBasis, 0, sizeof(gBasis));
	memset(ifBasis, 0, sizeof(ifBasis));
	RowFile.close();
	BasisFile.close();
	RowFile.open("wait.txt", ios::in | ios::out);
	BasisFile.open("nd.txt", ios::in | ios::out);
	//iToBasis.clear();

	ans.clear();
}

int readBasis() {          //读取消元子
	for (int i = 0; i < numBasis; i++) {
		if (BasisFile.eof()) return i - 1;
		string tmp;
		bool flag = false;
		int row = 0;
		getline(BasisFile, tmp);
		stringstream s(tmp);
		int pos;
		while (s >> pos) {
			if (!flag) {
				row = pos;
				flag = true;
				ifBasis[row] = 1;
			}
			int index = pos / 32;
			int offset = pos % 32;
			gBasis[row][index] = gBasis[row][index] | (1 << offset);
		}
		flag = false;
		row = 0;
	}
}

int readRowsFrom(int pos) {       //读取被消元行
	if (RowFile.is_open())
		RowFile.close();
	RowFile.open("wait.txt", ios::in | ios::out);
	memset(gRows, 0, sizeof(gRows));   //重置为0
	string line;
	for (int i = 0; i < pos; i++) {       //读取pos前的无关行
		getline(RowFile, line);
	}
	for (int i = pos; i < pos + maxrow; i++) {
		int tmp;
		getline(RowFile, line);
		if (line.empty()) return i;
		bool flag = false;
		stringstream s(line);
		while (s >> tmp) {
			int index = tmp / 32;
			int offset = tmp % 32;
			gRows[i - pos][index] = gRows[i - pos][index] | (1 << offset);
			flag = true;
		}
	}
	return -1;  //成功读取maxrow行

}

int findfirst(int row) {
	int first;
	for (int i = maxsize - 1; i >= 0; i--) {
		if (gRows[row][i] == 0)
			continue;
		else {
			int pos = i * 32;
			int offset = 0;
			for (int k = 31; k >= 0; k--) {
				if (gRows[row][i] & (1 << k))
				{
					offset = k;
					break;
				}
			}
			first = pos + offset;
			return first;
		}
	}
	return -1;
}

int _findfirst(int row) {
	int first;
	for (int i = maxsize - 1; i >= 0; i--) {
		if (answers[row][i] == 0)
			continue;
		else {
			int pos = i * 32;
			int offset = 0;
			for (int k = 31; k >= 0; k--) {
				if (answers[row][i] & (1 << k))
				{
					offset = k;
					break;
				}
			}
			first = pos + offset;
			return first;
		}
	}
	return -1;
}

void writeResult_MPI(ofstream& out) {
	for (int j = 0; j < num; j++) {
		for (int i = maxsize - 1; i >= 0; i--) {
			if (answers[j][i] == 0)
				continue;
			int pos = i * 32;
			for (int k = 31; k >= 0; k--) {
				if (answers[j][i] & (1 << k)) {
					out << k + pos << " ";
				}
			}
		}
		out << endl;
	}
}

void solve_MPI(int argc, char* argv[]) {
	int flag;
	double start_time = 0;
	double end_time = 0;
	MPI_Init(&argc, &argv);
	int total = 0;
	int rank = 0;
	int i = 0;
	int j = 0;
	int begin = 0, end = 0;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &total);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		flag = readRowsFrom(0);     //读取被消元行
		num = (flag == -1) ? maxrow : flag;
		begin = rank * num / total;
		end = (rank == total - 1) ? num : (rank + 1) * (num / total);
		for (i = 1; i < total; i++) {
			MPI_Send(&num, 1, MPI_INT, i, 0, MPI_COMM_WORLD);//0是被消元行行数
			int b = i * (num / total);
			int e = (i == total - 1) ? num : (i + 1) * (num / total);
			for (j = b; j < e; j++) {
				MPI_Send(&gRows[j][0], maxsize, MPI_INT, i, 1, MPI_COMM_WORLD);//1时被消元行数据
			}
		}

	}
	else {
		MPI_Recv(&num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		begin = rank * (num / total);
		end = (rank == total - 1) ? num : (rank + 1) * (num / total);
		for (i = begin; i < end; i++) {
			MPI_Recv(&gRows[i][0], maxsize, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime();
	for (i = begin; i < end; i++) {
		int first = findfirst(i);

		while (first != -1) {

			if (ifBasis[first] == 1) {
				for (j = 0; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ gBasis[first][j];

				}
				first = findfirst(i);
			}
			else {
				tmpAns.push_back(first);
				if (rank == 0) {
					for (j = 0; j < maxsize; j++) {
						gBasis[first][j] = gRows[i][j];
						answers[i][j] = gRows[i][j];
					}
					ifBasis[first] = 1;  //仅仅将0号进程消元到底
				}
				break;
			}
		}
		if (first == -1)
			tmpAns.push_back(-1);
	for (i = 0; i < rank; i++) {

		int b = i * (num / total);
		int e = b + num / total;
		for (j = b; j < e; j++) {
			MPI_Recv(&answers[j][0], maxsize, MPI_INT, i, 2, MPI_COMM_WORLD, &status);//接收来自进程i的消元结果，可能作为之后的消元子
			int first = _findfirst(j);
			firstToRow.insert(pair<int, int>(first, j));//记录下首项信息
		}
		for (j = begin; j < end; j++) {

			int first = tmpAns.at(j - begin);
			if (first == -1)
				continue;

			while ((firstToRow.find(first) != firstToRow.end() || ifBasis[first] == 1) && first != -1) {  //存在可消元项
				if (firstToRow.find(first) != firstToRow.end()) {  //消元结果有消元子
					int row = firstToRow.find(first)->second;
					for (int k = 0; k < maxsize; k++) {
						gRows[j][k] = gRows[j][k] ^ answers[row][k];
					}
					first = findfirst(i);
				}

				if (first == -1)
					break;
				if (ifBasis[first] == 1) {
					for (int k = 0; k < maxsize; k++) {
						gRows[j][k] = gRows[j][k] ^ gBasis[first][k];     //进行异或消元

					}
					first = findfirst(i);
				}
			}

		}
	}
	if (rank != 0) {
		for (i = begin; i < end; i++) {
			int first = findfirst(i);
			if (first == -1)
				continue;
			while ((firstToRow.find(first) != firstToRow.end() || ifBasis[first] == 1) && first != -1) {  //存在可消元项
				if (firstToRow.find(first) != firstToRow.end()) {
					int row = firstToRow.find(first)->second;
					for (int k = 0; k < maxsize; k++) {
						gRows[i][k] = gRows[i][k] ^ answers[row][k];
					}
					first = findfirst(i);
				}

				if (first == -1)
					break;
				if (ifBasis[first] == 1) {
					for (int k = 0; k < maxsize; k++) {
						gRows[i][k] = gRows[i][k] ^ gBasis[first][k];

					}
					first = findfirst(i);
				}
			}
			for (j = 0; j < maxsize; j++) {
				gBasis[first][j] = gRows[i][j];
				answers[i][j] = gRows[i][j];
			}
			ifBasis[first] = 1;


		}

	}
	for (i = rank + 1; i < total; i++) {
		for (j = begin; j < end; j++) {

			MPI_Send(&answers[j][0], maxsize, MPI_INT, i, 2, MPI_COMM_WORLD);
		}
	}

	if (rank == total - 1) {
		end_time = MPI_Wtime();
		cout << 1000 * (end_time - start_time) << "ms" << endl;
		writeResult_MPI(out_mpi);
		out_mpi.close();
	}
	MPI_Finalize();

}

int main(int argc, char* argv[]) {

	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

	readBasis();
	solve_MPI(argc, argv);
}
