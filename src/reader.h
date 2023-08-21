#ifndef READER_H
#define READER_H

#include <iostream>
#include<fstream>

#include <sys/time.h>

#include <unistd.h> // getcwd

using namespace std;



class Reader
{
    public:
        char CURRENT_DIR[500];
        char instanceG[50];
        char instanceKR[50];

        int num_vertices;/* numero de vertices/individuos do grafo*/
        int num_skills;/*numero de habilidades*/
        int num_teams; /* quantidade de total de equipes/projetos*/

        int **G; /*matriz |n|x|n| grafo*/
        int **K; /*matriz |n|x|f| de pessoa - habilidade*/
        int **R; /* matriz demanda |m|x|f|: quantidade de cada pessoa com habilidade k_a em cada equipe*/

        char* G_type;
    
    public:
        Reader(bool show, int vert, int graph_class, int class_type,const char* instance_type, char* graph_type="undirected"){
            
            this->G_type = graph_type;

            getcwd(CURRENT_DIR, 500);
            instanceReader(show,vert,graph_class,class_type,instance_type);
        }
        ~Reader(){
            // segmentation fault (core dump) error
            // for (int i = 0; i < num_vertices; i++){
            //     delete[] G[i];
            //     delete[] K[i];
            // }
            // for (int i = 0; i < num_skills; i++){
            //     delete[] R[i];
            // }
            num_vertices = 0; num_skills = 0; num_teams = 0;
        }
        // //--------------------------------------------read data------------------------------------------------------------------
        //load matrix G
        void load_Graph(char arq[], bool show, string instance_type){

            FILE *inst = fopen(arq,"r");
            if(!inst)
                printf("erro na leitura do arquivo!");

            int edge_weight;
            fscanf(inst,"%d",&num_vertices);
            G = (int**)(malloc(num_vertices*sizeof(int*)));

            for (int u=0;u<num_vertices;u++)
            {
                G[u] = (int*)malloc(num_vertices*sizeof(int));

                for (int v=0;v<num_vertices;v++)
                {
                    fscanf(inst,"%d",&edge_weight);
                    // G[u][v] = edge_weight;
                    if(u == v){
                        G[u][v] = 1;
                    }else{
                        if(edge_weight>0)
                            G[u][v] = 1;
                            // G[v][u] = 1;
                        else if(edge_weight<0)
                            G[u][v] = -1;
                            // G[v][u] = -1;
                        else if(edge_weight==0)
                            G[u][v] = 0;
                            // G[v][u] = 0;
                    }
                }

            }
            fclose(inst);


            if (instance_type == "random" && G_type == "undirected"){
                // cout << "random" << endl;
                for (int u = 0; u < num_vertices; u++){
                    for (int v = u+1; v < num_vertices; v++){
                        if (G[u][v] == 1 && G[v][u] != 1){
                            G[v][u] = 1;
                        }else if (G[v][u] == 1 && G[u][v] != 1){
                            G[u][v] = 1;
                        }else if (G[u][v] == -1 && G[v][u] != -1){
                            G[v][u] = -1;
                        }else if (G[v][u] == -1 && G[u][v] != -1){
                            G[u][v] = -1;
                        }
                    }
                }
            }


            if(show){
                cout << endl << "Numero de individuos: " << num_vertices << endl;
                cout << " Adjacency Matrix \n";
                for (int u = 0; u < num_vertices; u++){
                    for (int v = 0; v < num_vertices; v++)
                        cout <<  G[u][v] << " ";
                    cout<< "\n";
                }
            }
        }

        //load individuals skills
        void load_K(char arq[], bool show){

            FILE *inst = fopen(arq,"r");
            if(!inst)
                printf("erro na leitura do arquivo!");

            int value;
            fscanf(inst,"%d",&num_skills);
            K = (int**)(malloc(num_vertices*sizeof(int*)));
            for (int u=0;u<num_vertices;u++){
                K[u] = (int*)malloc(num_skills*sizeof(int));

                for (int s=0;s<num_skills;s++){
                    fscanf(inst,"%d",&value);
                    K[u][s] = value;
                }
            }
            if(show){
                printf("\nHabilidades = %d \n", num_skills);
                for (int u=0;u<num_vertices;u++){
                    for (int s=0;s<num_skills;s++)
                        printf("%d ", K[u][s]);
                    printf("\n");
                }
            }

            fclose(inst);
        }

        //load the skills necessary in projects/teams
        void load_R(char arq[], bool show){

            FILE *inst = fopen(arq,"r");
            if(!inst)
                printf("erro na leitura do arquivo!");

            int value;
            fscanf(inst,"%d",&num_teams);
            R = (int**)(malloc(num_teams*sizeof(int*)));

            for (int j=0;j<num_teams;j++){
                R[j] = (int*)malloc(num_skills*sizeof(int));

                for (int s=0;s<num_skills;s++){
                    fscanf(inst,"%d",&value);
                    R[j][s] = value;
                }
            }
            if(show){
                printf("\nEquipes = %d\n", num_teams);

                for (int j=0;j<num_teams;j++){
                    for (int s=0;s<num_skills;s++)
                        printf("%d ", R[j][s]);
                    printf("\n");
                }
            }

            fclose(inst);
        }
        
        void instanceReader(bool show, int vert, int graph_class, int class_type, string instance_type){
            char arq1[600]; char arq2[600]; char arq3[600];
            char arq_base[500];
            sprintf(arq_base, "%s/instances_COR", CURRENT_DIR);

            if (instance_type == "random"){
                sprintf(instanceG, "%dverticesS%d", vert, graph_class);
            }
            else if (instance_type == "bitcoinotc" || instance_type == "epinions"){
                if (G_type == "directed")
                    sprintf(instanceG, "%dvertices_%s_directed_S%d", vert, instance_type.c_str(), graph_class);
                else
                    sprintf(instanceG, "%dvertices_%s_S%d", vert, instance_type.c_str(), graph_class);

            }else{
                cout << "[INFO] not a valid instance type" << endl;
            }

            sprintf(arq1, "%s/%dVertices/%s.txt", arq_base, vert, instanceG);
            cout << "[INFO] Read graph: "<< instanceG << endl;
            cout<< "[INFO] Graph type: " << G_type << endl; 

            sprintf(instanceKR, "%dVertices/class1/%d", vert, class_type);
            sprintf(arq2, "%s/%s/K.txt", arq_base, instanceKR);
            sprintf(arq3, "%s/%s/R.txt", arq_base, instanceKR);
            cout << "[INFO] Instance class type: "<< instanceKR << endl;


            load_Graph(arq1, show, instance_type); load_K(arq2, show); load_R(arq3, show);

        }
    void show(){

        cout << endl << "Numero de individuos: " << num_vertices << endl;
                cout << " Adjacency Matrix \n";
                for (int u = 0; u < num_vertices; u++){
                    for (int v = 0; v < num_vertices; v++)
                        cout <<  G[u][v] << " ";
                    cout<< "\n";
                }

        printf("\nHabilidades = %d \n", num_skills);
                for (int u=0;u<num_vertices;u++){
                    for (int s=0;s<num_skills;s++)
                        printf("%d ", K[u][s]);
                    printf("\n");
                }

        printf("\nEquipes = %d\n", num_teams);

                for (int j=0;j<num_teams;j++){
                    for (int s=0;s<num_skills;s++)
                        printf("%d ", R[j][s]);
                    printf("\n");
                }
    }
};


#endif