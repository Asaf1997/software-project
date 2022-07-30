#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "spkmeans.h"


const int MAX_ITER_KMEANS = 300;
const int MAX_ITER_JACOBI = 100;


/* returns the distance between two vectors */
double distance(double vector1[], double vector2[], int length){
    double sum=0.0;
    int i;
    for(i=0; i < length; i++){
        sum += (vector1[i] - vector2[i])*(vector1[i] - vector2[i]);
    }
    sum = sqrt(sum);
    return sum;
}


void assertion_check(int b){
    if (b){
        printf("An Error Has Occurred");
        exit(1);
    }
}


void read_data(Graph * graph, char * file_path){
    double value;
    char c;
    int first_bool = 1, n = 1, d = 0, i, j;
    double ** vectors_array;
    FILE * f = fopen(file_path, "r");
    my_assert(f != NULL);

    while (fscanf(f, "%lf%c", &value, &c) == 2)
    {
        if (first_bool == 1) {
            d++;
        }

        if (c == '\n') {
            first_bool = 0;
            n++;
        }
    }

    rewind(f);

    vectors_array = make_mat(n, d);

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            fscanf(f, "%lf%c", &value, &c);
            vectors_array[i][j] = value;
        }
    }
    fclose(f);

    graph->vertices = vectors_array;
    graph->dim = d;
    graph->n = n;
    return;
}


double ** make_mat(int m, int n){
    double ** mat;
    int i,j;
    
    mat = calloc(m, sizeof(double*));
    assertion_check(mat == NULL);

    for (i = 0; i < m; i++){
        mat[i] = calloc(n, sizeof(double));
        assertion_check(mat[i] == NULL);

        for (j=0 ; j < n ; j++){
            mat[i][j] = 0.0;
        }
    }
    return mat;
}


double * make_vector(int n){
    double * vector = calloc(n, sizeof(double));
    assertion_check(vector == NULL);
    return vector;
}


void print_mat(double ** mat, int i, int j){
    int x,y;
    for (x = 0 ; x < i ; x++){
        for (y = 0 ; y < j ; y++){
            printf("%.4f", mat[x][y]);
            if (y != j-1){
                printf(",");
            }
        }
        printf("\n");
    }
}


void print_mat_transposed(double ** mat, int i, int j){
    int x,y;
    for (x = 0 ; x < j ; x++){
        for (y = 0 ; y < i ; y++){
            printf("%.4f", mat[y][x]);
            if (y != j-1){
                printf(",");
            }
        }
        printf("\n");
    }
}


void free_mat(double ** mat, int x){
    int i;
    for (i=0 ; i < x ; i++){
        free(mat[i]);
    }
    free(mat);
}


void wam(Graph * graph){
    double ** wam_mat;
    double ** vectors = graph->vertices;
    double weight;
    int i,j;
    int n = graph->n;

    wam_mat = make_mat(n, n);

    for (i = 0; i < n; i++){
        for (j=i+1 ; j < n ; j++){
            weight = exp( (-1) * distance(vectors[i], vectors[j], d) /2 );
            wam_mat[i][j] = weight;
            wam_mat[j][i] = weight;
        }
        wam_mat[i][i] = 0;
    }

    graph->vertices = wam_mat;
    return;
}


void ddg(Graph * graph){
    double ** weighted_adj_matrix = graph->wam_mat;
    double ** diagonal_matrix = make_mat(n, n);    
    double sum_weights = 0;
    int i;

    for(int i=0 ; i < n ; i++){
        sum_weights = 0;
        for(int j=0 ; j < n ; j++){
            sum_weights += weighted_adj_matrix[i][j];
        }
        diagonal_matrix[i][i] = sum_weights;
    }
    
    graph->ddg_mat = diagonal_matrix;
}


double ** sqrt_diagonal_matrix(double ** diagonal_matrix, int n){
    double ** sqrt_matrix = make_mat(n, n);

    for(int i=0 ; i < n ; i++){
        sqrt_matrix[i][i] = 1 / sqrt(diagonal_matrix[i][i]);
    }
    return sqrt_matrix;
}


static double ** lnorm(Graph * graph){
    double ** sqrt_D = sqrt_diagonal_matrix(ddg_mat, n);
    double ** lnorm_mat = make_mat(n, n);
    int i,j;

    if (sqrt_D == NULL || lnorm_mat == NULL){ return 1; }

    for (i=0 ; i < n ; i++){
        for(j=0 ; j < n ; j++){
            if (i==j){
                lnorm_mat[i][j] = 1.0;
            }
            lnorm_mat[i][j] = lnorm_mat[i][j] - wam_mat[i][j]*sqrt_D[i][i]*sqrt_D[j][j];
        }
    }
}


static double ** jacobi(double ** A, double ** eigenvectors, double * eigenvalues, int n){}


int main(int argc, char *argv[]){
    Graph graph = {0};
    Goal goal;
    int n, dim;
    char *input_name, *goal_str;
    double *eigenvalues;
    double **eigenvectors;

    /* checking correctness of the input */
    if(argc != 3 ){        
        printf("Invalid Input!");
        exit(1);
    }

    goal_str = argv[1];
    input_name = argv[2];

    if (goal_str != "wam" && goal_str != "ddg" && goal_str != "lnorm" && goal_str != "jacobi"){
        printf("Invalid Input!");
        exit(1);
    }

    goal = (int)goal_str[0];

    read_data(&graph, input_name);

    dim = graph.dim;
    n = graph.n;

    if (goal == e_jacobi){
        eigenvectors = make_mat(n ,n);
        eigenvalues = make_vector(n);
        jacobi(graph->vertices, eigenvectors, eigenvalues, n);

        print_mat(eigenvalues, 1 ,n);
        printf("\n");
        print_mat_transposed(eigenvectors, n, n);

        free(eigenvalues);
        free_mat(eigenvectors, n);
        free_mat(graph->vertices, n);
        free(graph);

        return 0;
    }

    wam(graph);
    if (goal == e_wam){
        print_mat(graph->wam_mat, n, n);
     
        free_mat(graph->vertices, n);
        free_mat(graph->wam_mat, n);
        free(graph);

        return 0;
    }

    ddg(graph);
    if (goal == e_ddg){
        print_mat(graph->ddg_mat, n, n);

        free_mat(graph->vertices, n);
        free_mat(graph->wam_mat, n);
        free_mat(graph->ddg_mat, n);
        free(graph);

        return 0;
    }

    lnorm(graph);

    print_mat(graph->lnorm_mat, n, n);

    free_mat(graph->vertices, n);
    free_mat(graph->wam_mat, n);
    free_mat(graph->ddg_mat, n);
    free_mat(graph->lnorm_mat, n);
    free(graph);

    return 0;
}