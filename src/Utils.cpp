#include "Utils.h"

using std::vector;

// comparacao entre doubles
int Utils::cmp(double x, double y, double tol) {
	return ( x <= y + tol ) ? ( x + tol < y ) ? -1 : 0 : 1;
}

// verifica se A eh um subconjunto de B, onde A armazena os seus elementos e B eh um vetor de booleanos indicando elementos pertencentes ao conjunto
bool Utils::isSubset(const vector< int > & A, const vector< int > & B) {
	for ( size_t i = 0; i < A.size(); ++i )	{
		if ( !B[ A[i] ] ) {
			return false;
		}
	}

	return true;
}

double Utils::roundPrecision(double value) {
	if(value >= -1.0*1e-6 && value <= 1e-6) {
		return 0.0;
	}
	else {
		return value;
	}
}

// Calcula a norma euclidiana de uma SEC representada por um vector v, de 0's e 1's, indicando os nós presentes ao conjunto S
double Utils::euclidianNorm(vector<int> &v, IloCplex * master_, IloNumArray &x_, vector<vector<int> > &edges_in_pool) {
	double sum = 0;
	int vsize = v.size();		// Número de vertices
	for(int i=0; i < vsize; i++) {
		for(int j=i+1; j < vsize; j++) {
			if(v[i] == 1 && v[j] == 1 && edges_in_pool[i][j] != -1) {		// Ambas as extremidades estão presentes no conj. S?
				sum += pow((double)x_[edges_in_pool[i][j]], 2);	// Adiciona à variável acumuladora sum o quadrado do valor de relaxação linear da aresta presente em S
			}
		}
	}
	return sqrt(sum);
}

double Utils::innerProduct(vector<int> v1, vector<int> v2) { // Assume-se que v1 e v2 tem o mesmo tamanho

	double ip = 0;
	int sz = v1.size();

	for(int i=0; i < sz; i++)
		ip += v1[i]*v2[i];

	return ip;
}


void Utils::bridgeUtil(int u, list<int> *adj, bool visited[], int disc[], int low[], int parent[], vector<vector<int> > &bridges) {

	// A static variable is used for simplicity, we can avoid use of static
	// variable by passing a pointer.
	static int time = 0;

	// Mark the current node as visited
	visited[u] = true;

	// Initialize discovery time and low value
	disc[u] = low[u] = ++time;

	// Go through all vertices aadjacent to this
	list<int>::iterator i;
	for (i = adj[u].begin(); i != adj[u].end(); ++i)
	{
		int v = *i;  // v is current adjacent of u

		// If v is not visited yet, then recur for it
		if (!visited[v])
		{
			parent[v] = u;
			bridgeUtil(v, adj, visited, disc, low, parent, bridges);

			// Check if the subtree rooted with v has a connection to
			// one of the ancestors of u
			low[u]  = min(low[u], low[v]);

			// If the lowest vertex reachable from subtree under v is
			// below u in DFS tree, then u-v is a bridge
			if (low[v] > disc[u]) {
				bridges[u][v] = bridges[v][u] = 1;
				//            	cout << u <<" " << v << endl;
			}
		}

		// Update low value of u for parent function calls.
		else if (v != parent[u])
			low[u]  = min(low[u], disc[v]);
	}
}

// DFS based function to find all bridges. It uses recursive function bridgeUtil()
void Utils::findBridges(vector<Edge> &tree, vector<vector<int> > &bridges) {

	// Mark all the vertices as not visited
	bool *visited = new bool[Instance::getInstance()->num_vertex];
	int *disc = new int[Instance::getInstance()->num_vertex];
	int *low = new int[Instance::getInstance()->num_vertex];
	int *parent = new int[Instance::getInstance()->num_vertex];

	// Initialize parent and visited arrays
	for (int i = 0; i < Instance::getInstance()->num_vertex; i++)
	{
		parent[i] = -1;
		visited[i] = false;
	}

	// Adjacency list
	list<int> *adj = new list<int>[Instance::getInstance()->num_vertex];
	for(unsigned int i=0; i<tree.size(); i++) {
		adj[tree[i].endpoint1].push_back(tree[i].endpoint2);
		adj[tree[i].endpoint2].push_back(tree[i].endpoint1);
	}

	// Call the recursive helper function to find Bridges
	// in DFS tree rooted with vertex 'i'
	for (int i = 0; i < Instance::getInstance()->num_vertex; i++)
		if (visited[i] == false)
			bridgeUtil(i, adj, visited, disc, low, parent, bridges);

	delete []adj;

}

void Utils::findHandles(vector<vector<int> > &adj_list, vector<vector<int> > &H) {

	char* checked = new char[Instance::getInstance()->num_vertex];		// o caracter w indica que o vértice ainda não foi visitado
	for(int i = 0; i < Instance::getInstance()->num_vertex; i++)
		checked[i] = 'w';

	for(int i = 0; i<Instance::getInstance()->num_vertex; i++) {
		if(checked[i] == 'w')
			H.push_back(bfsHandle(adj_list, checked, i));		// BFSHandle retorna um Handle e este é adicionado ao conjunto de Handles H
	}

}

vector<int> Utils::bfsHandle(vector<vector<int> > &adj_list, char checked[], int source) {

	vector<int> h(Instance::getInstance()->num_vertex, 0); 	// Conjunto Handle - H[j] = 1 indica que o vértice j
	// pertence ao handle, e H[j] = 0 caso contrário
	queue<int> Q;
	checked[source] = 'b';
	h[source] = 1;
	Q.push(source);
	while(Q.empty() == false) {
		int u = Q.front();
		Q.pop();
		for(unsigned int i=0; i < adj_list[u].size(); i++) {
			int v = adj_list[u][i];
			if(checked[v] == 'w') {
				checked[v] = 'b';
				h[v] = 1;
				Q.push(v);
			}
		}
	}

	return h;
}

void Utils::bfs(int source, char checked[], vector<vector<int> > &lista_adj) {

	queue<int> Q;
	checked[source] = 'b';
	Q.push(source);

	while(Q.empty() == false) {
		int u = Q.front();
		Q.pop();
		for(int i=0; i< (int)lista_adj[u].size(); i++) {
			int v = lista_adj[u][i];
			if(checked[v] == 'w') {
				checked[v] = 'b';
				Q.push(v);
			}
		}
	}
}

void Utils::checa_solucao_bfs (vector<Edge> v) {


	vector<vector<int> > lista_adj(Instance::getInstance()->num_vertex);
	char* checked = new char[Instance::getInstance()->num_vertex];
	for(int i = 0; i<Instance::getInstance()->num_vertex; i++)
		checked[i] = 'w';

	for(int i=0; i<(int)v.size(); i++) {
		lista_adj[v[i].endpoint1].push_back(v[i].endpoint2);
		lista_adj[v[i].endpoint2].push_back(v[i].endpoint1);
	}

	vector<int> s;
	int cont = 0;
	for(int i = 0; i<Instance::getInstance()->num_vertex; i++) {
		if(checked[i] == 'w') {
			cont++;
			Utils::bfs(i, checked, lista_adj);
		}
	}
	cout << "# comp. conexas " << cont << endl;
	cin.get();
}

