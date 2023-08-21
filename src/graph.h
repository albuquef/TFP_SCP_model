#ifndef GRAPH_H
#define GRAPH_H



#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>


#include <iostream>
#include <limits>
#include <algorithm>

#include "matching.h"

using namespace std;

#define MAXV 1000
#define INF 1000000 // std::numeric_limits<int>::max();

class Graph{
    public:
        int **G; // matrix graph = v,e
        int **w; // matrix graph weight :  default = [1]*num_vertices
        int num_vertices;// num of vertices |V|
        char* G_type; // directed or undirected
        bool connected;
        int **Graph_SPP; // matrix graph weighted by th Shortest Positive Path
        float **Even_Graph; // signed graph to even graph 
        int cont_pos_edges;
        int **Edmonds_Graph_st; // shortest even path to min weighted matching
        float **Weighted_Edmond_Graph_st; // shortest even path to min weighted matching
    public:
        Graph(int **G,int num_vertices,char* G_type);
        ~Graph();
        void printGraph();
        void traverse(int u, bool visited[]);
        bool isConnected();
        void genereateGraph_weighted(char* method);
        void generateEven_Graph();
        double runEdmonds_Graph_st(int s, int t);

};

class EdgeNode{
    public:
        int key;
        int weight;
        int sign; // +1(default) or -1  ps: is possible add the concept of +weight and -weight?
        EdgeNode *next;
        EdgeNode(int key, int weight,int sign=1){
            this->key = key;
            this->weight = weight;
            this->next = NULL;
            this->sign = sign;
        }
};

class GraphAdj{
    // bool directed;
    // int num_vertices;
    public:
        bool directed;
        int num_vertices;
        EdgeNode *edges[MAXV];
        GraphAdj(bool directed,int num_vertices=MAXV){
            this->directed = directed;
            this->num_vertices=num_vertices;
            for(int i = 0; i < (num_vertices); i ++){
                this->edges[i] = NULL;
            }
        }
        ~GraphAdj(){//to do
        };
        void insert_edge(int x, int y, int weight, bool directed){
            if(x >= 0 && x < (num_vertices) && y >= 0 && y < (num_vertices)){
                EdgeNode *edge = new EdgeNode(y, weight);
                edge->next = this->edges[x];
                this->edges[x] = edge;
                if(!directed){
                    insert_edge(y, x, weight, true);
                }
            }
        }
        void insert_edge_sign(int x, int y, int weight,int sign, bool directed){
            if(x >= 0 && x < (num_vertices) && y >= 0 && y < (num_vertices)){
                EdgeNode *edge = new EdgeNode(y, weight,sign);
                edge->next = this->edges[x];
                this->edges[x] = edge;
                if(!directed){
                    insert_edge_sign(y, x, weight,sign, true);
                }
            }
        }
        void print(){
            for(int v = 0; v < (num_vertices); v ++){
                if(this->edges[v] != NULL){
                    cout << "Vertex " << v << " has neighbors: " << endl;
                    EdgeNode *curr = this->edges[v];
                    while(curr != NULL){
                        cout << curr->key << endl;
                        curr = curr->next;
                    }
                }
            }
        }
        void init_vars(bool discovered[], int distance[], int parent[]){
            for(int i = 0; i < (num_vertices); i ++){
                discovered[i] = false;
                distance[i] = INF;
                parent[i] = -1;
            }
        }
        void dijkstra_shortest_path(GraphAdj *g, int parent[], int distance[], int start){

            bool discovered[num_vertices];
            EdgeNode *curr;
            int v_curr;
            int v_neighbor;
            int weight;
            int smallest_dist;

            init_vars(discovered, distance, parent);

            distance[start] = 0;
            v_curr = start;

            while(discovered[v_curr] == false){

                discovered[v_curr] = true;
                curr = g->edges[v_curr];

                while(curr != NULL){ //look the distance of one edge/arc from the v_curr (neighbor) 

                    v_neighbor = curr->key;
                    weight = curr->weight;

                    if((distance[v_curr] + weight) < distance[v_neighbor]){
                        distance[v_neighbor] = distance[v_curr] + weight;
                        parent[v_neighbor] = v_curr;
                    }
                    curr = curr->next;
                }

                //set the next current vertex to the vertex with the smallest distance
                smallest_dist = std::numeric_limits<int>::max();
                for(int i = 0; i < (num_vertices); i ++){
                    if(!discovered[i] && (distance[i] < smallest_dist)){
                        v_curr = i;
                        smallest_dist = distance[i];
                    }
                }
            }
        }
        
        void init_vars_sign(bool discovered_pos[],bool discovered_neg[], int parent_pos[], int parent_neg[], int distance_pos[], int distance_neg[]){
            for(int i = 0; i < (num_vertices); i ++){
                discovered_pos[i] = false;
                discovered_neg[i] = false;
                distance_pos[i] = INF;
                distance_neg[i] = INF;
                parent_pos[i] = -1;
                parent_pos[i] = -1;
            }
        }
        void dijkstraSIGN_shortest_path(GraphAdj *g, int parent_pos[], int parent_neg[], int distance_pos[], int distance_neg[], int start){

            bool discovered_pos[num_vertices]; 
            bool discovered_neg[num_vertices];
            EdgeNode *curr;
            int v_curr;
            int v_neighbor;
            int weight;
            int sign;
            int smallest_dist_pos;
            int smallest_dist_neg;

            init_vars_sign(discovered_pos, discovered_neg, parent_pos, parent_neg, distance_pos, distance_neg);

            distance_pos[start] = 0;
            distance_neg[start] = INF;
            v_curr = start;

            while(discovered_pos[v_curr] == false){

                discovered_pos[v_curr] = true;
                discovered_neg[v_curr] = true;
                curr = g->edges[v_curr];
            
                while(curr != NULL){ //look the distance of one edge/arc from the v_curr (neighbor) 

                    v_neighbor = curr->key;
                    weight = curr->weight;
                    sign = curr->sign;
                    
                    if (sign == 1 && v_neighbor != start){
                        
                        if((distance_pos[v_curr] + weight) < distance_pos[v_neighbor]){
                            distance_pos[v_neighbor] = distance_pos[v_curr] + weight;
                            parent_pos[v_neighbor] = v_curr;
                        } else if((distance_neg[v_curr] + weight) < distance_neg[v_neighbor]){
                            distance_neg[v_neighbor] = distance_neg[v_curr] + weight;
                            parent_neg[v_neighbor] = v_curr;
                        }
                    }
                    else if (sign == -1 && v_neighbor != start){
                        
                        if((distance_pos[v_curr] + weight) < distance_neg[v_neighbor]){
                            distance_neg[v_neighbor] = distance_pos[v_curr] + weight;
                            parent_neg[v_neighbor] = v_curr;
                        }
                        if((distance_neg[v_curr] + weight) < distance_pos[v_neighbor]){
                            distance_pos[v_neighbor] = distance_neg[v_curr] + weight;
                            parent_pos[v_neighbor] = v_curr;
                        }

                        
                    }
                    curr = curr->next;
                }

                                   //set the next current vertex to the vertex with the smallest distance
                    int smallest_dist = INF;
                    for(int i = 0; i < (num_vertices); i ++){
                        if(!discovered_pos[i] && (distance_pos[i] < smallest_dist)){
                            v_curr = i;
                            smallest_dist = distance_pos[i];
                        }
                        if(!discovered_neg[i] && (distance_neg[i] < smallest_dist)){
                            v_curr = i;
                            smallest_dist = distance_neg[i];
                        }
                        
                    }     
            }
        }
        void print_shortest_path(int v, int parent[]){
            if(v >= 0 && v < (num_vertices) && parent[v] != -1){
                cout << parent[v] << " ";
                print_shortest_path(parent[v], parent);
            }
        }
        void print_distances(int start, int distance[]){
            for(int i = 0; i < (num_vertices); i ++){
                if(distance[i] < INF){
                    cout << "Shortest distance from " << start << " to " << i << " is: " << distance[i] << endl;
                }
            }
        /*
         teste dijkstra:
                GraphAdj *g = new GraphAdj(false,6);
                int parent[MAXV + 1];
                int distance[MAXV + 1];
                int start = 1;

                g->insert_edge(1, 2, 4, false);
                g->insert_edge(1, 3, 1, false);
                g->insert_edge(3, 2, 1, false);
                g->insert_edge(3, 4, 5, false);
                g->insert_edge(2, 4, 3, false);
                g->insert_edge(2, 5, 1, false);
                g->insert_edge(4, 5, 2, false);

                g->dijkstra_shortest_path(g, parent, distance, start);
                //print shortest path from vertex 1 to 5
                g->print_shortest_path(5, parent);
                g->print_distances(start, distance);
        */



}

};

#endif // GRAPH_H