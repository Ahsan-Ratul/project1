#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define MAX_GENES 100
#define MAX_CONDITIONS 10

float data[MAX_GENES][MAX_CONDITIONS];       
float graph[MAX_GENES][MAX_GENES];           
int rows = 0, cols = 0;                    
int parent[MAX_GENES];                       
int visited[MAX_GENES];                       


int cluster_id[MAX_GENES];

void read_data_from_input() {
    printf("Enter number of genes (rows): ");
    scanf("%d", &rows);
    printf("Enter number of conditions (columns): ");
    scanf("%d", &cols);

    printf("Enter gene expression values (each row separated by new line):\n");
    for (int i = 0; i < rows; i++) {
        printf("Gene %d: ", i + 1);
        for (int j = 0; j < cols; j++) {
            scanf("%f", &data[i][j]);
        }
    }
}

float euclidean_distance(float *a, float *b, int length) {
    float sum = 0;
    for (int i = 0; i < length; i++) {
        float diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}


void build_distance_matrix() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            graph[i][j] = euclidean_distance(data[i], data[j], cols);
        }
    }
}


void print_distance_matrix() {
    printf("\nDistance Matrix (Euclidean):\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            printf("%.2f ", graph[i][j]);
        }
        printf("\n");
    }
}


void save_distance_matrix_to_file(const char *filename) {
    FILE *fout = fopen(filename, "w");
    if (fout == NULL) {
        printf("❌ Error: Could not open output file.\n");
        return;
    }

    fprintf(fout, "Distance Matrix (Euclidean):\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            fprintf(fout, "%.2f ", graph[i][j]);
        }
        fprintf(fout, "\n");
    }

    fclose(fout);
    printf("✅ Distance matrix saved to '%s'\n", filename);
}


void prim_mst() {
    float key[MAX_GENES];
    int i, count, u, v;

    for (i = 0; i < rows; i++) {
        key[i] = FLT_MAX;
        visited[i] = 0;
        parent[i] = -1;
    }

    key[0] = 0;

    for (count = 0; count < rows - 1; count++) {
        float min = FLT_MAX;
        u = -1;

        for (i = 0; i < rows; i++) {
            if (!visited[i] && key[i] < min) {
                min = key[i];
                u = i;
            }
        }

        visited[u] = 1;

        for (v = 0; v < rows; v++) {
            if (graph[u][v] != 0 && !visited[v] && graph[u][v] < key[v]) {
                parent[v] = u;
                key[v] = graph[u][v];
            }
        }
    }

    printf("\nMinimum Spanning Tree (MST):\n");
    float total_cost = 0;
    for (i = 1; i < rows; i++) {
        printf("Gene %d -- Gene %d   Distance: %.2f\n", parent[i] + 1, i + 1, graph[i][parent[i]]);
        total_cost += graph[i][parent[i]];
    }
    printf("Total MST Cost: %.2f\n", total_cost);
}


void save_mst_to_file(const char *filename) {
    FILE *fout = fopen(filename, "w");
    if (fout == NULL) {
        printf("❌ Error: Could not open MST output file.\n");
        return;
    }

    fprintf(fout, "EdgeList (MST):\n");
    for (int i = 1; i < rows; i++) {
        fprintf(fout, "%d %d %.2f\n", parent[i] + 1, i + 1, graph[i][parent[i]]);
    }
    fclose(fout);
    printf("✅ MST edges saved to '%s'\n", filename);
}


int find_cluster(int x) {
    if (cluster_id[x] != x)
        cluster_id[x] = find_cluster(cluster_id[x]);
    return cluster_id[x];
}

void union_clusters(int a, int b) {
    int rootA = find_cluster(a);
    int rootB = find_cluster(b);
    if (rootA != rootB) {
        cluster_id[rootB] = rootA;
    }
}


void cluster_genes(float threshold) {
    
    for (int i = 0; i < rows; i++) {
        cluster_id[i] = i;
    }

    
    for (int i = 1; i < rows; i++) {
        float dist = graph[i][parent[i]];
        if (dist <= threshold) {
            union_clusters(i, parent[i]);
        }
    }

    
    int cluster_num = 0;
    int cluster_map[MAX_GENES] = {0}; 

    printf("\nClusters formed with threshold %.2f:\n", threshold);
    for (int i = 0; i < rows; i++) {
        int root = find_cluster(i);
        if (cluster_map[root] == 0) {
            cluster_num++;
            cluster_map[root] = cluster_num;
        }
        printf("Gene %d -> Cluster %d\n", i + 1, cluster_map[root]);
    }

    printf("Total clusters formed: %d\n", cluster_num);
}

int main() {
    read_data_from_input();                 
    build_distance_matrix();               
    print_distance_matrix();               
    save_distance_matrix_to_file("distance_matrix.txt"); 
    prim_mst();
    save_mst_to_file("mst_edges.txt");

    float threshold;
    printf("\nEnter distance threshold to form clusters (e.g., 2.5): ");
    scanf("%f", &threshold);
    cluster_genes(threshold);

    return 0;
}







