#ifndef MATCHING_H
#define MATCHING_H


#include <list>
#include <vector>

#include "graph.h"
using namespace std;

// class BinaryHeap
#define EPSILON 0.000001
#define INFINITO 1000000000.0
#define GREATER(A, B) ((A) - (B) > EPSILON)
#define LESS(A, B) ((B) - (A) > EPSILON)
#define EQUAL(A, B) (fabs((A) - (B)) < EPSILON)
#define GREATER_EQUAL(A, B) (GREATER((A),(B)) || EQUAL((A),(B)))
#define LESS_EQUAL(A, B) (LESS((A),(B)) || EQUAL((A),(B)))
#define MIN(A, B) (LESS((A),(B)) ? (A) : (B))
#define MAX(A, B) (LESS((A),(B)) ? (B) : (A))


// class Matching 
#define EVEN 2
#define ODD 1
#define UNLABELED 0


class Graph_Matching
{
public:
	//n is the number of vertices
	//edges is a list of pairs representing the edges (default = empty list)
	Graph_Matching(int n, const list< pair<int, int> > & edges = list< pair<int, int> >()):
                    n(n),
                    m(edges.size()),
                    adjMat(n, vector<bool>(n, false)),
                    adjList(n),
                    edges(),
                    edgeIndex(n, vector<int>(n, -1)){

                    for(list< pair<int, int> >::const_iterator it = edges.begin(); it != edges.end(); it++)
                    {
                        int u = (*it).first;
                        int v = (*it).second;

                        AddEdge(u, v);
                    }
    }

	//Default constructor creates an empty Graph_Matching
	Graph_Matching(): n(0), m(0) {};

	//Returns the number of vertices
	int GetNumVertices() const { return n; };
	//Returns the number of edges
	int GetNumEdges() const { return m; };

	//Given the edge's index, returns its endpoints as a pair
	pair<int, int> GetEdge(int e) const{
    	if(e > (int)edges.size())
	    	throw "Error: edge does not exist";

	    return edges[e];
    }
	//Given the endpoints, returns the index
	int GetEdgeIndex(int u, int v) const{
        if( u > n or
            v > n )
            throw "Error: vertex does not exist";

        if(edgeIndex[u][v] == -1)
            throw "Error: edge does not exist";

        return edgeIndex[u][v];
    }

	//Adds a new vertex to the Graph_Matching
	void AddVertex(){

        for(int i = 0; i < n; i++)
        {
            adjMat[i].push_back(false);
            edgeIndex[i].push_back(-1);
        }
        n++;
        adjMat.push_back( vector<bool>(n, false) );
        edgeIndex.push_back( vector<int>(n, -1) );
        adjList.push_back( list<int>() );
    }
    
	//Adds a new edge to the Graph_Matching
	void AddEdge(int u, int v){
        
        if( u > n or
            v > n )
            throw "Error: vertex does not exist";

        if(adjMat[u][v]) return;

        adjMat[u][v] = adjMat[v][u] = true;
        adjList[u].push_back(v);
        adjList[v].push_back(u);

        edges.push_back(pair<int, int>(u, v));
        edgeIndex[u][v] = edgeIndex[v][u] = m++;

    }

	//Returns the adjacency list of a vertex
	const list<int> & AdjList(int v) const
    {
        if(v > n)
            throw "Error: vertex does not exist";

        return adjList[v];
    }

	//Returns the Graph_Matching's adjacency matrix
	const vector< vector<bool> > & AdjMat() const{
    	return adjMat;
    }
private:
	//Number of vertices
	int n;
	//Number of edges
	int m;

	//Adjacency matrix
	vector< vector<bool> > adjMat;

	//Adjacency lists
	vector< list<int> > adjList;

	//Array of edges
	vector< pair<int, int> > edges;

	//Indices of the edges
	vector< vector<int> > edgeIndex;
};


class BinaryHeap
{
public:
	BinaryHeap(): satellite(1), size(0) {};

	//Inserts (key k, satellite s) in the heap
	void Insert(double k, int s){

        //Ajust the structures to fit new data
        if(s >= (int)pos.size())
        {
            pos.resize(s+1, -1);
            key.resize(s+1);
            //Recall that position 0 of satellite is unused
            satellite.resize(s+2);
        }
        //If satellite is already in the heap
        else if(pos[s] != -1)
        {
            throw "Error: satellite already in heap";
        }

        int i;
        for(i = ++size; i/2 > 0 && GREATER(key[satellite[i/2]], k); i /= 2)
        {
            satellite[i] = satellite[i/2];
            pos[satellite[i]] = i;
        }
        satellite[i] = s;
        pos[s] = i;
        key[s] = k;
    }
	//Deletes the element with minimum key and returns its satellite information
	int DeleteMin(){

        if(size == 0)
            throw "Error: empty heap";

        int min = satellite[1];
        int slast = satellite[size--];


        int child;
        int i;
        for(i = 1, child = 2; child  <= size; i = child, child *= 2)
        {
            if(child < size && GREATER(key[satellite[child]], key[satellite[child+1]]))
                child++;

            if(GREATER(key[slast], key[satellite[child]]))
            {
                satellite[i] = satellite[child];
                pos[satellite[child]] = i;
            }
            else
                break;
        }
        satellite[i] = slast;
        pos[slast] = i;

        pos[min] = -1;

        return min;
    }
	//Changes the key of the element with satellite s
	void ChangeKey(double k, int s){
        Remove(s);
        Insert(k, s);
    }
	//Removes the element with satellite s
	void Remove(int s){
        int i;
        for(i = pos[s]; i/2 > 0; i /= 2)
        {
            satellite[i] = satellite[i/2];
            pos[satellite[i]] = i;
        }
        satellite[1] = s;
        pos[s] = 1;

        DeleteMin();
    }
	//Returns the number of elements in the heap
	int Size(){
        return size;
    }
	//Resets the structure
	void Clear()
    {
        key.clear();
        pos.clear();
        satellite.clear();
    }

private:
	vector<double> key;//Given the satellite, this is its key
	vector<int> pos;//Given the satellite, this is its position in the heap
	vector<int> satellite;//This is the heap!

	//Number of elements in the heap
	int size;
};


class Matching
{
public:
	//Parametric constructor receives a graph instance
	Matching(const Graph_Matching & G);
    // Matching(int **G);

	//Solves the minimum cost perfect matching problem
	//Receives the a vector whose position i has the cost of the edge with index i
	//If the graph doest not have a perfect matching, a const char * exception will be raised
	//Returns a pair
	//the first element of the pair is a list of the indices of the edges in the matching
	//the second is the cost of the matching
	pair< list<int>, double > SolveMinimumCostPerfectMatching(const vector<double> & cost);

	//Solves the maximum cardinality matching problem
	//Returns a list with the indices of the edges in the matching
	list<int> SolveMaximumMatching();

private:
	//Grows an alternating forest
	void Grow();
	//Expands a blossom u
	//If expandBlocked is true, the blossom will be expanded even if it is blocked
	void Expand(int u, bool expandBlocked);
	//Augments the matching using the path from u to v in the alternating forest
	void Augment(int u, int v);
	//Resets the alternating forest
	void Reset();
	//Creates a blossom where the tip is the first common vertex in the paths from u and v in the hungarian forest
	int Blossom(int u, int v);
	void UpdateDualCosts();
	//Resets all data structures 
	void Clear();
	void DestroyBlossom(int t);
	//Uses an heuristic algorithm to find the maximum matching of the graph
	void Heuristic();
	//Modifies the costs of the graph so the all edges have positive costs
	void PositiveCosts();
	list<int> RetrieveMatching();

	int GetFreeBlossomIndex();
	void AddFreeBlossomIndex(int i);
	void ClearBlossomIndices();

	//An edge might be blocked due to the dual costs
	bool IsEdgeBlocked(int u, int v);
	bool IsEdgeBlocked(int e);
	//Returns true if u and v are adjacent in G and not blocked
	bool IsAdjacent(int u, int v);

	const Graph_Matching & G;

	list<int> free;//List of free blossom indices

	vector<int> outer;//outer[v] gives the index of the outermost blossom that contains v, outer[v] = v if v is not contained in any blossom
	vector< list<int> > deep;//deep[v] is a list of all the original vertices contained inside v, deep[v] = v if v is an original vertex
	vector< list<int> > shallow;//shallow[v] is a list of the vertices immediately contained inside v, shallow[v] is empty is the default
	vector<int> tip;//tip[v] is the tip of blossom v
	vector<bool> active;//true if a blossom is being used

	vector<int> type;//Even, odd, neither (2, 1, 0)
	vector<int> forest;//forest[v] gives the father of v in the alternating forest
	vector<int> root;//root[v] gives the root of v in the alternating forest 

	vector<bool> blocked;//A blossom can be blocked due to dual costs, this means that it behaves as if it were an original vertex and cannot be expanded
	vector<double> dual;//dual multipliers associated to the blossoms, if dual[v] > 0, the blossom is blocked and full
	vector<double> slack;//slack associated to each edge, if slack[e] > 0, the edge cannot be used
	vector<int> mate;//mate[v] gives the mate of v

	int m, n;

	bool perfect;

	list<int> forestList;
	vector<int> visited;
};


#endif // MATCHING_H