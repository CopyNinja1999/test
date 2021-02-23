
/*
Going deeper:1.When the cache miss happened
2.how large the size is the algorithm become ineffictive
3.why using deque is faster
4.the performance in HPC Cluster against the performance in desktop
5.

*/

#include <iostream>
#include <fstream>
#include <stack>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cmath>
#include <deque>
#include <queue>
#include <string>
#include <sstream>
//#include "Vertice.h"
//#include "Graph.h";
using namespace std;
using std::vector;
using std::pair;
using std::stack;//neccessary
using std::queue;
struct vertice {
	size_t postOrder = 0;
	size_t preOrder = 0;
	size_t ccNum = 0;
	size_t index = 0;
	bool visited = false;
};
size_t cc = 1;//count isolated components
size_t preCounter,postCounter;
vertice* v;//record the status
vector<vector<int> > adj;//record the neighbours
//void printV(vertice v) {
//	size_t i = v.index;
//	size_t n = adj[i].size();
//	if (n == 0) {
//		puts("");
//		return;
//	}
//	else {
//		for (size_t j = 0; j < n; j++) {
//			cout << adj[i][j] << " ";
//		}
//		puts("");
//	}
//}
bool equal(vector<int> &temp,size_t x) {
	//if one veritce adjecancy list temp owns the same element x, this function will return true
	size_t n = temp.size();
	bool equal = false;
	if (n == 0) {
		return false;
	}
	else {
		for (size_t c = 0; c < n; c++)
		{
			if (c >= n || c < 0) { 
				cout << "vetcor out of range" << endl; 
			break; }
			if (temp[c] == x)
			{
				equal = true;
				return equal;
			}

		}
		return equal;
	}
}
size_t Num_of_edges(vector<vector<int>> &adj) {
	size_t m=0;
	size_t n = adj.size();
	//in the adjacancy list, the number of edges is the sum of the edges that each vertice owns
	for (size_t i = 0; i < n;i++) {
		m += adj[i].size();
	}


	return m;
}
void Print_adj_list() {
	int n = adj.size();
	for (int i = 0; i < n; i++) {
		int m = adj[i].size();
		cout << i << ": ";
		for (int j = 0; j < m; j++) { cout << adj[i][j] << " ,"; }
		cout << endl;

	}

}
void Print_adj_matrix() {
	int n = adj.size();
	for (int i = 0; i < n;i++) {//outer loop, print i rows
		for (int j = 0; j < n; j++) {//inner loop, print dots
			if (equal(adj[i], j)) //which means it has the element j
			{
				cout << "1" << " ";
			}
			else
				cout << "0" << " ";
		}
		cout << endl;
	
	}
}
size_t random(int r, int b) {
	size_t random;
	random=(rand() % (b-r+1)) +r;//represents [r,b]

//	cout << random << endl;
	return random;
}
bool isExplored() {
	size_t n = adj.size();
	bool eFlag =true;
	for (size_t i = 0; i < n; i++)
	{
		if (!v[i].visited) {
			eFlag = false;
			cout << "Traversal incomplete"<<endl;
			break;
		}
	}
	
	return eFlag;
}
//check if the vector temp has the same element as x

void generateD(size_t n, size_t m)
{
	
	cout << "Directed graph generated:" << endl;
	cout << "n=" << n << " m=" << m << endl;
	adj = vector<vector<int> >(n, vector<int>());
	v = new vertice[n];
	for (size_t i = 0; i < n; i++) {
		v[i].index = i;
	}
	for (size_t i = 0; i < m; i++)//push m lines in template vector
	{
		size_t j, k;
		j = random(1, n); //[1,n]
		while (true) {
			if (adj[j - 1].size() == n - 1)
			{
				j = random(1, n);
			}
			else
				break;
		}
		k = random(1, n);
		while (true) {
			if (j == k) {
				k = random(1, n);

			}
			else break;
		}
		if (!equal(adj[j - 1], k - 1))  // if elements are not repeatd or are empty, return false, push element to the back
		{
			adj[j - 1].push_back(k - 1);
		}
		else {
			//make sure there are exactly m edges in the graph, in case the vector has the same elements.
			while (true)
			{
				k = random(1, n);
				if (!equal(adj[j - 1], k - 1)&&k!=j)//
				{
					adj[j - 1].push_back(k - 1);
					break;
				}
			};

		}
	}
	//if (n < 25) {
	//	for (size_t i = 0; i < n; i++)
	// {
	//	cout << i << ": ";
	//	printV(v[i]);

	// }
	//}

}
void generateU(size_t n, size_t m) {
	cout << "Undirected graph generated:" << endl;
	cout << "n=" << n << " m=" << m << endl;
	adj = vector<vector<int> >(n, vector<int>());
	v = new vertice[n];
	for (size_t i = 0; i < n; i++) {
		v[i].index = i;
	}
	for (size_t i = 0; i < m; i++)//push m lines in template vector
	{
		size_t j, k;
		j = random(1, n); //[1,n]
		//for the smaller example, in case all edges are on the same vertice
		while (true) 
		{
			if (adj[j - 1].size() == n - 1)
			{
				j = random(1, n);
			}
			else
				break;
		}
		k = random(1, n);
		while (true) {
			if (j == k) {
				k = random(1, n);

			}
			else break;
		}
		if (!equal(adj[j - 1], k - 1) && !equal(adj[k - 1], j - 1))  // if elements are not repeatd or are empty, return false, push element to the back
		{
			adj[j - 1].push_back(k - 1);
			adj[k - 1].push_back(j - 1);
		}
		else {
			//make sure there are exactly m edges in the graph, in case the vector has the same elements.
			while (true)
			{
				k = random(1, n);
				if (!equal(adj[j - 1], k - 1) && !equal(adj[k - 1], j - 1)&&k != j)
				{
					adj[j - 1].push_back(k - 1);
					adj[k - 1].push_back(j - 1);
					break;
				}
			};

		}
	}
	//if (n < 25) {
	//	for (size_t i = 0; i < n; i++)
	//	{
	//		cout << i << ": ";
	//		printV(v[i]);

	//	}
	//}

}
void generateBi(size_t n, size_t m) {
	cout <<"bipartite graph generated"<< endl;
	cout << "n=" << n << " m=" << m << endl;
	//distribute n1 vertices for the first block and m1 edges towards the other block
	size_t n1 = random(1,n);
	size_t m1 = random(1,m);
	adj = vector<vector<int> >(n, vector<int>());
	v = new vertice[n];
	for (size_t i = 0; i < n; i++) {
		v[i].index = i;
	}
	//assign edges from the first block
	for (size_t i = 0; i < m1; i++)//push m1 lines in template vector
	{
		size_t j, k;
		j = random(1, n1); //[1,n1]
		k = random(n1+1,n);//[n1+1,n]
		if (!equal(adj[j - 1], k - 1) )  // if elements are not repeatd or are empty, return false, push element to the back
		{
			adj[j - 1].push_back(k - 1);
		}
		else {
			//make sure there are exactly m edges in the graph, in case the vector has the same elements.
			while (true)
			{
				k = random(n1+1,n);
				if (!equal(adj[j - 1], k - 1))
				{
					adj[j - 1].push_back(k - 1);
					break;
				}
			}


		}
	}
//and the second block
	for (size_t i = 0; i < m-m1; i++)//push m1 lines in template vector
	{
		size_t j, k;
		j = random(1, n1); //[1,n1]
		k = random(n1 + 1, n);//[n1+1,n]
		if (!equal(adj[k - 1], j - 1))  // if elements are not repeatd or are empty, return false, push element to the back
		{
			adj[k - 1].push_back(j - 1);
		}
		else {
			//make sure there are exactly m edges in the graph, in case the vector has the same elements.
			while (true)
			{
				j = random(1, n1);
				if (!equal(adj[k - 1], j - 1))
				{
					adj[k - 1].push_back(j - 1);
					break;
				}
			}


		}
	}
}
void userInterface() 
{
	size_t a;
	size_t n=0, m=0;
	cout << "Welcome to algorithms test, please select the type of the graph:" << "\n 1.small 2.large\n" << endl;
	cin >> a;
	switch (a) 
	{
	case 1:
		cout << "You selected small.\n" << endl;
	    n = random(5, 20);
		break;
	case 2:
		cout << "You selected small.\n" << endl;
		n = random(1000, 100000);
		break;
	default:
		cout << "Input illegal.\n" << endl;
		break;
	}	
	cout << " 1.sparse 2.dense\n" << endl;
	cin >> a;
	switch (a)
	{
	case 1:
		cout << "You selected sparse.\n"<<endl;
		m = n+2;
		break;
	case 2:
		cout << "You selected dense.\n" << endl;
		m= (size_t)(pow(n,2)-n-2)/2;//undirected
		break;
	default:
		cout << "Input illegal.\n" << endl;
		break;
	}

	cout << " 1.directed 2.undirected\n" << endl;
	size_t b;
	cin >> b;
	switch (b) {
	case 1:
		
		cout << "You selected directed.\n" << endl;
		generateD(n, m);
		
			break;
	case 2:
		cout << "You selected undirected.\n" << endl;
		generateU(n, m);
		break;
	default:
		cout << "Input illegal.\n" << endl;
		break;

		}
	cout << "Print the graph?(Recommended when small graph is selected)" << "\n 1.Yes 2.No\n" << endl;
	cin >> a;
	switch (a)
	{
	case 1:
		//for (size_t i = 0; i < n; i++)
		//{
		//	cout << i << ": ";
		//	printV(v[i]);

		//}
		break;
	case 2:
		break;
	default:
		cout << "Input illegal.\n" << endl;
		break;
	}

	};



void explore(int i)
{
	if (!v[i].visited)
	{
		v[i].visited = true;
		v[i].ccNum = cc;
		int n = adj[i].size();
		if (!adj[i].empty())
		{
			for (int j = 0; j < n; j++)
			{int w = adj[i][j];//acutal index in the vector
				explore(w);
			}
			v[i].postOrder = postCounter++;
		}
		else { 
			v[i].postOrder = postCounter++;
			return; }//if the vertice is sink, return the upper level 
	}


	else {
		
		return;} //if visited, return immediately
}
void dfsRecursive() //the naive implementation of directed dfs
{
	postCounter = 0;
	int n = adj.size();
	for (int w = 0; w < n; w++) {
		if (!v[w].visited) {
			explore(w);
			cc++;
		}
	}

}	
void dfsStack() {//using naive c++ stl stack
	size_t n = adj.size();
	for (size_t i = 0; i < n; i++)
	{
		v[i].visited = false;
	}
	stack<vertice> S;
	for (size_t k = 0; k < n; k++) {
		if (v[k].visited) { continue; }
		else
			v[k].visited = true;//start point v[0]
		S.push(v[k]);
		while (!S.empty())
		{
			vertice u = S.top();
			S.pop();//these two steps represent the operation pop
			size_t m = adj[u.index].size();
			int  index = -1;
			for (size_t j = 0; j < m; j++)
			{
				size_t i = adj[u.index][j];
				if (!v[i].visited)
				{
					index = i;
					break;
				}
			}
			if (index != -1)
			{
				S.push(u);
				v[index].visited = true;
				S.push(v[index]);
			}
		}
	}


	S.~stack<vertice>();//release memory
}
void dfsDeque() {
	size_t n = adj.size();
	for (size_t i = 0; i < n; i++)
	{
		v[i].visited = false;
	}
	deque<vertice> D;
	for (size_t k = 0; k < n; k++) {
		if (v[k].visited) { continue; }
		else
			v[k].visited = true;//start point v[0]
		D.push_front(v[k]);
		while (!D.empty())
		{
			vertice u = D.front();
			D.pop_front();
			size_t m = adj[u.index].size();
			int index = -1;
			for (size_t j = 0; j < m; j++)
			{
				size_t i = adj[u.index][j];
				if (!v[i].visited)
				{
					index = i;
					break;
				}
			}
			if (index != -1)
			{
				D.push_front(u);
				v[index].visited = true;
				D.push_front(v[index]);
			}
		}
	}
	D.~deque<vertice>();//release memory

}
void showOrder() {
	size_t n = adj.size();
	for (size_t i = 0; i < n; i++)
	{
		cout << i << ": " << v[i].postOrder << endl;

		
	}
}
void bfsQueue() 
{
	size_t n = adj.size();
	for (size_t i = 0; i < n; i++)
	{
		v[i].visited = false;
	}
	queue<vertice> Q;
	for (size_t k = 0; k < n; k++) {
		if (v[k].visited) { 
			continue; 
		}
		else
			v[k].visited = true;
		Q.push(v[k]);
		while (!Q.empty())
		{
			vertice u = Q.front();
			Q.pop();//these two steps represent the operation dequeue()
			size_t m = adj[u.index].size();
			for (size_t j = 0; j < m; j++)
			{
				size_t i = adj[u.index][j];
				if (!v[i].visited)
				{
					v[i].visited = true;
					Q.push(v[i]);
				}
			}

		}
	}
	Q.~queue<vertice>();
}
void bfsDeque() {
	size_t n = adj.size();
	for (size_t i = 0; i < n; i++)
	{
		v[i].visited = false;
	}
	deque<vertice> D;
	for (size_t k = 0; k < n; k++) {
		if (v[k].visited) { continue; }
		else
			v[k].visited = true;
		D.push_back(v[k]);
		while (!D.empty())
		{
			vertice u = D.front();
			D.pop_front();//these two steps represent the operation dequeue()
			size_t m = adj[u.index].size();
			for (size_t j = 0; j < m; j++)
			{
				size_t i = adj[u.index][j];
				if (!v[i].visited)
				{
					v[i].visited = true;
					D.push_back(v[i]);
				}
			}

		}
	}
	D.~deque<vertice>();
}
void test_directed() {
	ofstream exp_data("dir_test.csv");// For the convenience of analyzing
	size_t n, m;

	for (size_t i = 1; i < 100; i++) {
		n = 100 * i;
		//m = n + 2;
		m = (size_t)(pow(n, 2) - n - 2) / 2;// m ->n^2
		generateD(n,m);
		clock_t t1, t2, t3;
		t1 = clock();//get current time
		dfsStack();
		t2 = clock();
		double dif1, dif2;
		dif1 = difftime(t2, t1);
		dfsDeque();
		t3 = clock();
		dif2 = difftime(t3, t2);
		exp_data << dif1 <<","<< dif2 << "\n";
	}
	exp_data.close();
}
int Max(size_t a, size_t b, size_t c)
{
	size_t max;
	if (a >= b)
	{
		if (a >= c) {
			max = a;
		}
		else
			max = c;
	}
	else if (b >= c) { max = b; }
	else max = c;
	return max;
}
void compareDFS() {
	clock_t t1,t2,t3;
	t1 = clock();//get current time
	dfsStack();
	t2 = clock();
	double dif1,dif2;
	dif1 = difftime(t2, t1);
	dfsDeque();
	t3 = clock();
	dif2 = difftime(t3, t2);
	cout << "Naive stack takes " <<dif1<<" ms" <<endl;
	cout << "Deque takes " << dif2 << " ms" << endl;
}
void compareBFS() {
	clock_t t1, t2, t3;
	t1 = clock();//get current time
	bfsQueue();
	t2 = clock();
	double dif1, dif2;
	dif1 = difftime(t2, t1);
	bfsDeque();
	t3 = clock();
	dif2 = difftime(t3, t2);
	cout << "Naive queue takes " << dif1 << " ms" << endl;
	cout << "Deque takes " << dif2 << " ms" << endl;
}
void Input_data(const string& filename) {
	std::fstream in(filename.c_str());
	cout << "reading file:" << filename << endl;
	string s;
	size_t n = 0, m = 0;
	string data1, data2;
	while (true)
	{
		std::getline(in, s);
		istringstream is(s);
		is >> data1 >> data2;
		int d1 = stoi(data1);
		int d2 = stoi(data2);
		n = Max(n, d2, d1);
		m += 1;
		if (in.eof()) { break; }
	}
	//this block will count the number of lines and calculate the maximun number appeared in the file, which are the parameters n, m(vertice, edge)
	in.seekg(0, ios::beg);
	n += 1;
	adj = vector<vector<int> >(n, vector<int>());
	for (size_t i = 0; i < m; i++)
	{
		int x, y;
		std::getline(in, s);
		istringstream is(s);
		is >> data1 >> data2;
		x = stoi(data1);
		y = stoi(data2);
		adj[x].push_back(y);
	}
	in.close();
	//this block will assign data into the vertice template in terms of the adjancancy list
}


void Output_data() {
	ofstream data("out.txt");
	size_t m = Num_of_edges(adj);
	for (size_t i = 0; i < m; i++)
	{
		size_t n = adj[i].size();
		for (size_t j = 0; j < n; j++) 
		{		data <<i<<" "<< adj[i][j] << "\n";}

	}
	data.close();
}


int main(){

	//Input_data("out.txt");

//srand((int)time(NULL));  // generate random seeds, use and only use once, ensure that every execuation is random
//generateBi(1480000,1480000);	//size_t range 18446744073709551615(2^64-1), int range 2,147,483,647
////Print_adj_matrix();
//	Print_adj_list();	
	//compareDFS();
//
//	isExplored();
//	Output_data();
	test_directed();
	return 0;
}