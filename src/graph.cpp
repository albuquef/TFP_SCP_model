#include "graph.h"



Graph::Graph(int **G,int num_vertices,char* G_type="undirected"){

   getcwd(CURRENT_DIR, 500);

   this->G = G;
   this->G_type = G_type;
   this->num_vertices = num_vertices;
   if(isConnected()){
      cout<< "[INFO] Graph is connected" << endl;
      this->connected=true;
   }else{
      cout<< "[INFO] Graph is not connected" << endl;
      this->connected=false;
   }

   // genereateGraph_weighted_paths("DijkstraComp");
   // genereateGraph_weighted_paths("MinMatching");


   // generateEdmonds_Graph_st(0,2);
   //  GraphAdj *g = new GraphAdj(false,6);


   // int parent_pos[MAXV + 1];
   // int distance_pos[MAXV + 1];
   // int parent_neg[MAXV + 1];
   // int distance_neg[MAXV + 1];

   // g->insert_edge_sign(1, 2, 4, -1, false);
   // g->insert_edge_sign(1, 3, 1, 1, false);
   // g->insert_edge_sign(3, 2, 1, -1, false);
   // g->insert_edge_sign(3, 4, 5, -1, false);
   // g->insert_edge_sign(2, 4, 3, -1, false);
   // g->insert_edge_sign(2, 5, 1, -1, false);
   // g->insert_edge_sign(4, 5, 2, 1, false);


   // g->dijkstraSIGN_shortest_path(g, parent_pos, parent_neg, distance_pos,distance_neg, start);
   //  //print shortest path from vertex 1 to 5
   // int final = 4;
   // cout << "Path reverse: "<< final << " ";
   // // g->print_shortest_path(final, parent_pos);
   // for(int i =0; i<5;i++){
   //    cout << parent_pos[i] << " ";
   // }
   // cout << endl; 
   // for(int i =0; i<5;i++){
   //    cout << parent_neg[i] << " ";
   // }
   // cout << endl; 
   // // g->print_distances(start, distance_pos);


}


Graph::~Graph(){

}

void Graph::printGraph(){
      cout << endl << "Numero de individuos: " << num_vertices << endl;
      cout << " Adjacency Matrix  \n" ;
      for (int u = 0; u < num_vertices; u++){
         for (int v = 0; v < num_vertices; v++)
            cout << G[u][v] << " ";
         cout<< "\n";
      }
}

void Graph::traverse(int u, bool visited[]){
   visited[u] = true; //mark v as visited
   for(int v = 0; v<num_vertices; v++){
      if(G[u][v]){
         if(!visited[v])
         traverse(v, visited);
      }
   }
}
bool Graph::isConnected(){
   bool *vis = new bool[num_vertices];
   //for all vertex u as start point, check whether all num_verticess are visible or not
   for(int u; u < num_vertices; u++){
      for(int i = 0; i<num_vertices; i++)
      vis[i] = false; //initialize as no num_vertices is visited
      traverse(u, vis);
      for(int i = 0; i<num_vertices; i++){
         if(!vis[i]) //if there is a num_vertices, not visited by traversal, graph is not connected
         return false;
      }
   }
   return true;
}


void Graph::generateEven_Graph(){

   cont_pos_edges = 0;
   for (int u = 0; u < num_vertices; u++){
      for (int v = u+1; v < num_vertices; v++)
            if(G[u][v] > 0)
               cont_pos_edges += 1;
      }

   // Even_Graph = malloc((num_vertices+cont_pos_edges)*sizeof(float*));
   // for (int u=0;u<num_vertices;u++)
   // {
   //    Even_Graph[u] = malloc((num_vertices+cont_pos_edges)*sizeof(float));
   // }

   Even_Graph = new float * [(num_vertices+cont_pos_edges)];
   for (int u = 0; u < (num_vertices+cont_pos_edges); u++) {
      Even_Graph[u] = new float[(num_vertices+cont_pos_edges)];
   }


   float a = 0.0;
   for (int u = 0; u < (num_vertices+cont_pos_edges); u++)
      for (int v = 0; v < (num_vertices+cont_pos_edges); v++){
         Even_Graph[u][v]= a;
   }

   int cont = num_vertices;
   for (int u = 0; u < num_vertices; u++)
      for (int v = u+1; v < num_vertices; v++){

            // cout << u+1 << " , " << v+1 << endl;
            if(G[u][v] < 0){
               Even_Graph[u][v] = static_cast<float>(abs(G[u][v]));
               Even_Graph[v][u] = static_cast<float>(abs(G[u][v]));
            }
            if (G[u][v] > 0){
               Even_Graph[u][cont] = static_cast<float>(G[u][v])/2;
               Even_Graph[cont][u] = static_cast<float>(G[u][v])/2;

               Even_Graph[cont][v] = static_cast<float>(G[u][v])/2;
               Even_Graph[v][cont] = static_cast<float>(G[u][v])/2;

               cont += 1;
            }
            if (cont > (num_vertices+cont_pos_edges)){
               cout << "error: cont > num_vertices + cont_pos_edges" << endl;
            }
         }
}
double Graph::runEdmonds_Graph_st(int s, int t){

   // Edmonds Graph st = H-s, H-t, barra(E)
   int num_vertices_EG = num_vertices+cont_pos_edges;
   // int num_vertices_edmonds = 2*(num_vertices_EG);
   int num_vertices_edmonds = 2*(num_vertices_EG) + 2; // +2 auxiliary vertices
   // cout << num_vertices_edmonds<< endl;
   

   // int **Edmonds_Graph_st; // shortest even path to min weighted matching
   // float **Weighted_Edmond_Graph_st; // shortest even path to min weighted matching


   // alocate memory
   int ** Edmonds_Graph_st = new int * [num_vertices_edmonds];
   for (int u = 0; u < num_vertices_edmonds; u++) {
      Edmonds_Graph_st[u] = new int[num_vertices_edmonds];
   }

   // alocate memory
   float ** Weighted_Edmond_Graph_st = new float * [num_vertices_edmonds];
   for (int u = 0; u < num_vertices_edmonds; u++) {
      Weighted_Edmond_Graph_st[u] = new float[num_vertices_edmonds];
   }

   // incialize
   for(int u=0;u<num_vertices_edmonds;u++)
      for(int v=0;v<num_vertices_edmonds;v++){
         Edmonds_Graph_st[u][v] = 0.0;
         Weighted_Edmond_Graph_st[u][v] = 0.0;
      }

   for (int u = 0; u < num_vertices_EG; u++)
      for (int v = u+1; v < num_vertices_EG; v++){

            if (u+num_vertices_EG > num_vertices_edmonds || v+num_vertices_EG > num_vertices_edmonds){
               cout << "error: position > num_vertices_edmonds" << endl;
            }
            if(u != t && v!= t && abs(Even_Graph[u][v])){
               // cout << "u',v' = " << u+1 <<"," << v+1 <<endl;
               Edmonds_Graph_st[u][v] = 1;
               Edmonds_Graph_st[v][u] = 1;

               Weighted_Edmond_Graph_st[u][v] = abs(Even_Graph[u][v]);
               Weighted_Edmond_Graph_st[v][u] = abs(Even_Graph[u][v]);

            }
            if (u != s && v!= s && abs(Even_Graph[u][v])){
               // cout << "u'',v'' = " << u+num_vertices_EG+1 <<"," << v+num_vertices_EG+1 <<endl;
            
               Edmonds_Graph_st[u+num_vertices_EG][v+num_vertices_EG] = 1;
               Edmonds_Graph_st[v+num_vertices_EG][u+num_vertices_EG] = 1;

               Weighted_Edmond_Graph_st[u+num_vertices_EG][v+num_vertices_EG] = abs(Even_Graph[u][v]);
               Weighted_Edmond_Graph_st[v+num_vertices_EG][u+num_vertices_EG] = abs(Even_Graph[u][v]);
            
            }
            if(u==s || v ==s){
               Edmonds_Graph_st[u+num_vertices_EG][v+num_vertices_EG] = 0; Edmonds_Graph_st[v+num_vertices_EG][u+num_vertices_EG] = 0;

               Weighted_Edmond_Graph_st[u+num_vertices_EG][v+num_vertices_EG] = 10000000;
               Weighted_Edmond_Graph_st[v+num_vertices_EG][u+num_vertices_EG] = 10000000;
            }
            if(u==t || v ==t){
               Edmonds_Graph_st[u][v] = 0; Edmonds_Graph_st[v][u] = 0;

               Weighted_Edmond_Graph_st[u][v] = 10000000;
               Weighted_Edmond_Graph_st[v][u] = 10000000;
            }
         }


      // edges \bar{E} u' and u'' for all V\{s',t''}
      for(int u=0; u<num_vertices_EG;u++){
         if( u!=s  && u!=t){
            // cout << "u', u'' = " << u+1 << "," << u+num_vertices_EG+1 << endl;
            Edmonds_Graph_st[u][u+num_vertices_EG] = 1;
            Weighted_Edmond_Graph_st[u][u+num_vertices_EG] = 0; //0.0000001;

            Edmonds_Graph_st[u+num_vertices_EG][u] = 1;
            Weighted_Edmond_Graph_st[u+num_vertices_EG][u] = 0;//0.0000001;
         }else{
            Edmonds_Graph_st[u][u+num_vertices_EG] = 0;
            Weighted_Edmond_Graph_st[u][u+num_vertices_EG] = 10000000; //0.0000001;

            Edmonds_Graph_st[u+num_vertices_EG][u] = 0;
            Weighted_Edmond_Graph_st[u+num_vertices_EG][u] = 10000000;//0.0000001;
         }
      }

      // gambiarra para sempre ter emparelhamento para os vértices s'' e t' que são criados sem necessidade
      // ex 4verticesS2 |V| = 12 + 2     s,t = 0,3 => 6 - (w = 0) - 12 - w(inf) - 0 e  8 - (w = 0) - 13 - w(inf) - 3    
      // s + num_vertices_EG - (w = 0) - num_vertices_edmonds-2 - w(inf) - s | t + num_vertices_EG - (w = 0) - num_vertices_edmonds-1 - w(inf) - s
      // outras arestas para os vertices auxiliares (num_vertices_edmonds-2) e (num_vertices_edmonds-1) não existem
      
      int v_aux_1 = num_vertices_edmonds-2; int v2 = s+num_vertices_EG; 
      Edmonds_Graph_st[v2][v_aux_1] = 1; Edmonds_Graph_st[v_aux_1][v2] = 1; 
      Weighted_Edmond_Graph_st[v2][v_aux_1] = 0; Weighted_Edmond_Graph_st[v_aux_1][v2] = 0; //0.0000001;

      Edmonds_Graph_st[v_aux_1][s] = 1; Edmonds_Graph_st[s][v_aux_1] = 1;
      Weighted_Edmond_Graph_st[v_aux_1][s] = 10000000;
      Weighted_Edmond_Graph_st[s][v_aux_1] = 10000000;


      int v_aux_2= num_vertices_edmonds-1; int v4 = t+num_vertices_EG;
      Edmonds_Graph_st[v4][v_aux_2] = 1; Edmonds_Graph_st[v_aux_2][v4] = 1;
      Edmonds_Graph_st[v4][v_aux_2] = 10000000; 
      Edmonds_Graph_st[v_aux_2][v4] = 10000000;

      Edmonds_Graph_st[v_aux_2][t] = 1; Edmonds_Graph_st[t][v_aux_2] = 1;
      Weighted_Edmond_Graph_st[v_aux_2][t] = 0; Weighted_Edmond_Graph_st[t][v_aux_2] = 0; //0.0000001;



      // cout << endl << "Numero de individuos: " << num_vertices_edmonds << endl;
      // cout << " Adjacency Matrix Edmonds Graph: (s,t) = (" << s << "," << t << ") \n" ;
      // for (int u = 0; u < num_vertices_edmonds; u++){
      //    for (int v = 0; v < num_vertices_edmonds; v++)
      //       cout << Edmonds_Graph_st[u][v] << " ";
      //    cout<< "\n";
      // }

      int num_edges_edmonds = 0;
      for (int u = 0; u < num_vertices_edmonds; u++){
         for (int v = u+1; v < num_vertices_edmonds; v++)
            if(Edmonds_Graph_st[u][v] > 0)
               num_edges_edmonds+=1;
         // cout<< "\n";
      }

   // cout << "------------" << endl;


   //    cout << endl << "Numero de individuos: " << num_vertices_edmonds << endl;
   //    cout << " Adjacency Matrix Edmonds Graph: (s,t) = (" << s << "," << t << ") \n" ;
   //    for (int u = 0; u < num_vertices_edmonds; u++){
   //       for (int v = 0; v < num_vertices_edmonds; v++)
   //          cout << Weighted_Edmond_Graph_st[u][v] << " ";
   //       cout<< "\n";
   //    }

      Graph_Matching G_match(num_vertices_edmonds);
      vector<double> cost(num_edges_edmonds+2);
      for (int u = 0; u < num_vertices_edmonds; u++){
         for (int v = u+1; v < num_vertices_edmonds; v++){
            if(Edmonds_Graph_st[u][v]>0){
               // cout << "add edge: " << u+1 << "," << v+1 <<endl;
               G_match.AddEdge(u, v);
               cost[G_match.GetEdgeIndex(u, v)] = Weighted_Edmond_Graph_st[u][v];
               // cost[G_match.GetEdgeIndex(u, v)] = 1;
            }
         }
      }


   // cout << "num vertices: " << num_vertices_edmonds << endl; 
   // cout << "num vertices edmonds: " << G_match.GetNumVertices() << endl;
	// cout << "\n";
   // cout << "num edges: " << num_edges_edmonds << endl; 
   // cout << "num edges edmonds: " << G_match.GetNumEdges() << endl;
   
   //Create a Matching instance passing the graph
	Matching M(G_match);


	// //Pass the costs to solve the problem
	pair< list<int>, double > solution = M.SolveMinimumCostPerfectMatching(cost);

	list<int> matching = solution.first;
	double obj = solution.second;

	// cout << "Optimal matching cost: " << obj << endl;
	// cout << "Edges in the matching:" << endl;
	// for(list<int>::iterator it = matching.begin(); it != matching.end(); it++)
	// {
	// 	pair<int, int> e = G_match.GetEdge( *it );

	// 	cout << e.first+1 << " " << e.second+1 << endl;
	// }


   // free memory matrix and matching

   for(int i = 0; i < num_vertices_edmonds; ++i)
  {
     delete Edmonds_Graph_st[i];
     delete Weighted_Edmond_Graph_st[i];
  }

  delete[] Edmonds_Graph_st;
  delete[] Weighted_Edmond_Graph_st;



   return obj;
}


void Graph::genereateGraph_weighted_paths(const char* method)
{
   cout << "[INFO] Creating graph weighted by SCPs" << endl;
   cout << "[INFO] Method: " << method << endl;
   
   Graph_SPP = (int**)(malloc(num_vertices*sizeof(int*)));
   for (int u=0;u<num_vertices;u++)
   {
      Graph_SPP[u] = (int*)malloc(num_vertices*sizeof(int));
   }

   if (method=="DijkstraComp" || strcmp(method,"DijkstraComp") == 0){

      cout << "[WARNING] Signed Graph must be balanced " << endl; // todo: check if graph is not balanced

      GraphAdj *g;
      if (G_type == "directed"){
         g = new GraphAdj(true,num_vertices);
      }else{
         g = new GraphAdj(false,num_vertices);
      }

      // matrix to adj list edges
      for(int u = 0;u<num_vertices;u++){
         for(int v = 0;v<num_vertices;v++){
            if(G[u][v] > 0){
               g->insert_edge_sign(u, v, G[u][v], 1, g->directed); // g->directed = true or false
            }else if(G[u][v] < 0){
               g->insert_edge_sign(u, v, abs(G[u][v]), -1, g->directed);
            }
         }
      }

      for(int u=0;u<num_vertices;u++){
         int parent_pos[num_vertices];
         int distance_pos[num_vertices];
         int parent_neg[num_vertices];
         int distance_neg[num_vertices];

         g->dijkstraSIGN_shortest_path(g, parent_pos, parent_neg, distance_pos,distance_neg, u);

         for(int v=0;v<num_vertices;v++){
            Graph_SPP[u][v] = distance_pos[v];
         }
      }
   } 
   else if(method == "MinMatching" || strcmp(method,"MinMatching") == 0){

      cout << "[WARNING] Signed Graph must be undirected " << endl;  // todo: check if graph is not directed
      // create Even Graph Gp
      generateEven_Graph();

      for(int u=0;u<num_vertices;u++){
         Graph_SPP[u][u] = 0;
      }


      // create and run Edmonds Graph with Min Weighted Perfect Matching
      double value = -1;
      for(int u=0;u<num_vertices;u++)
         for(int v=u+1; v<num_vertices;v++){
            try
            {
               value = runEdmonds_Graph_st(u,v);
               Graph_SPP[u][v] = int(value);
               Graph_SPP[v][u] = int(value);
            }
            catch(const char * msg) // not exist perfect matching => not exist positive path
            {  
               Graph_SPP[u][v] = 0;
               Graph_SPP[v][u] = 0;
               // cout << msg << endl;
               // return 1;
            }

         }



   }


}

void Graph::saveResults_SPP(string instance_type,int gclass,const char* method, double timeTotal){

   char instanceG[50];

   if (instance_type == "random"){
         sprintf(instanceG, "%dverticesS%d", num_vertices, gclass);
         if (strcmp(G_type,"directed") == 0)
            sprintf(instanceG, "%dverticesS%d_directed", num_vertices, gclass);
   }
   else if (instance_type == "bitcoinotc" || instance_type == "epinions"){
         if (strcmp(G_type,"directed") == 0){
            sprintf(instanceG, "%dvertices_%s_directed_S%d", num_vertices, instance_type.c_str(), gclass);
         }
         else{
            sprintf(instanceG, "%dvertices_%s_S%d", num_vertices, instance_type.c_str(), gclass);
         }
   }else{
         cout << "[INFO] not a valid instance type" << endl;
   }



   char arq[1000];
    // sprintf(arq, "%s/results/%d_Vertices_2022-06-27_cycle_mtz.ods",CURRENT_DIR, rd->num_vertices);
    if (strcmp(G_type,"directed") == 0)
        sprintf(arq, "%s/results/result_%d_Vertices_directed_SPP_%s.ods",CURRENT_DIR, num_vertices,method);
    else
        sprintf(arq, "%s/results/result_%d_Vertices_SPP_%s.ods",CURRENT_DIR, num_vertices,method);
    

    cout << arq << endl;

    ofstream outputTable;
    outputTable.open(arq,ios:: app);
    if(outputTable.is_open()){

        outputTable << instanceG << ";"; // grafo instancia
        outputTable << num_vertices << ";";   // numero de vertices
        outputTable << timeTotal <<  ";"; // tempo execucao tfp (cplex)
        outputTable << " \n ";


    }

    outputTable.close();











}