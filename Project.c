// gene_expression_clustering_programiz.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_GENES 100
#define MAX_CONDITIONS 10

float data[MAX_GENES][MAX_CONDITIONS];
int rows = 0, cols = 0;

// --- Read data from standard input (works in online compilers like Programiz) ---
void read_data_from_input() {
    printf("Enter number of genes (rows): ");
    scanf("%d", &rows);
    printf("Enter number of conditions (columns): ");
    scanf("%d", &cols);

    printf("Enter gene expression values (each row separated by new line):\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            scanf("%f", &data[i][j]);
        }
    }
}

// --- Distance Calculation (Euclidean) ---
float euclidean_distance(float *a, float *b, int length) {
    float sum = 0;
    for (int i = 0; i < length; i++) {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(sum);
}

void print_distance_matrix() {
    printf("\nDistance Matrix (Euclidean):\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            float d = euclidean_distance(data[i], data[j], cols);
            printf("%.2f ", d);
        }
        printf("\n");
    }
}

int main() {
    read_data_from_input();
    print_distance_matrix();
    return 0;
}

/*
Sample Input (copy-paste into Programiz input box):
3
3
1.2 2.3 3.4
2.1 3.2 4.3
1.1 2.0 3.1
*/


