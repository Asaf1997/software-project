#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "spkmeans.h"


const int MAX_ITER_KMEANS = 300;
const int MAX_ITER_JACOBI = 100;


int main(int argc, char *argv[]){
    Graph * graph = malloc(sizeof(Graph*));
    Goal goal;
    int n;
    char *input_name, *goal_str;
    double *eigenvalues;
    double **eigenvectors;

    /* checking correctness of the input */
    if(argc != 3 ){
        printf("Invalid Input!\n");
        exit(1);
    }
    
    goal_str = argv[1];
    input_name = argv[2];
    

    if (strcmp(goal_str, "wam") &&  strcmp(goal_str, "ddg") && strcmp(goal_str, "lnorm") && strcmp(goal_str, "jacobi")){
        printf("Invalid Input!\n");
        exit(1);
    }

    goal = (int)goal_str[0];

    /* insert information from file to graph */
    read_data(graph, input_name);

    n = graph->n;

    /* operate by goal */
    if (goal == e_jacobi){
        eigenvectors = make_mat_identity(n ,n);
        eigenvalues = make_vector(n);
        jacobi(graph->vertices, n, eigenvectors, eigenvalues);

        print_array(eigenvalues, n);
        print_mat(eigenvectors, n, n);

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
    int first_bool = 1, n = 0, d = 0, i, j;
    double ** vectors_array;
    FILE * f = fopen(file_path, "r");
    assertion_check(f == NULL);

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


double ** make_mat_identity(int m, int n){
    double ** mat;
    int i,j;
    
    mat = calloc(m, sizeof(double*));
    assertion_check(mat == NULL);

    for (i = 0; i < m; i++){
        mat[i] = calloc(n, sizeof(double));
        assertion_check(mat[i] == NULL);

        for (j=0 ; j < n ; j++){
            if (i == j) { mat[i][j] = 1.0; } else { mat[i][j] = 0.0; }
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


void print_array(double * array, int n){
    int y;
    for (y = 0 ; y < n ; y++){
            printf("%.4f", array[y]);
            if (y != n-1){
                printf(",");
            }
    }
    printf("\n");
}


void free_mat(double ** mat, int x){
    int i;
    if (mat != NULL){
        for (i=0 ; i < x ; i++){
        free(mat[i]);
        }
    }
    free(mat);
}


void wam(Graph * graph){
    double ** wam_mat;
    double ** vectors = graph->vertices;
    double weight;
    int i,j;
    int n = graph->n;
    int d = graph->dim;

    wam_mat = make_mat(n, n);

    for (i = 0; i < n; i++){
        for (j=i+1 ; j < n ; j++){
            weight = exp( (-1) * (distance(vectors[i], vectors[j], d) /2) );
            wam_mat[i][j] = weight;
            wam_mat[j][i] = weight;
        }
        wam_mat[i][i] = 0;
    }

    graph->wam_mat = wam_mat;
    return;
}


void copy_mat(double ** old, double ** copy, int rows, int culs){
    int i,j;
    copy = make_mat(rows, culs);
    for(i=0 ; i<rows ; i++){
        for(j=0 ; j<culs ; j++){
            copy[i][j] = old[i][j];
        }
    }
}


void ddg(Graph * graph){
    double ** weighted_adj_matrix = graph->wam_mat;
    double ** diagonal_matrix;    
    double sum_weights = 0;
    int i,j;
    int n = graph->n;

    diagonal_matrix = make_mat(n, n); 

    for(i=0 ; i < n ; i++){
        sum_weights = 0;
        for(j=0 ; j < n ; j++){
            sum_weights += weighted_adj_matrix[i][j];
        }
        diagonal_matrix[i][i] = sum_weights;
    }
    
    graph->ddg_mat = diagonal_matrix;
}


double ** sqrt_diagonal_matrix(double ** diagonal_matrix, int n){
    double ** sqrt_matrix = make_mat(n, n);
    int i;
    for(i=0 ; i < n ; i++){
        sqrt_matrix[i][i] = 1 / sqrt(diagonal_matrix[i][i]);
    }
    return sqrt_matrix;
}


void lnorm(Graph * graph){
    double ** sqrt_D = sqrt_diagonal_matrix(graph->ddg_mat, graph->n);
    double ** lnorm_mat;
    double ** wam_mat = graph->wam_mat;
    int i,j;
    int n = graph->n;
    lnorm_mat = make_mat(n,n);

    for (i=0 ; i < n ; i++){
        for(j=0 ; j < n ; j++){
            if (i==j){
                lnorm_mat[i][j] = 1.0;
            }
            lnorm_mat[i][j] = lnorm_mat[i][j] - wam_mat[i][j]*sqrt_D[i][i]*sqrt_D[j][j];
        }
    }

    graph->lnorm_mat = lnorm_mat;
    free_mat(sqrt_D, n);
    return;
}


void find_largest_element_pivot(double ** matrix, int n, int * piv){
    double max_val = 0.0;
    int i,j;
    for (i = 0; i < n; i++){
        for (j = i + 1 ; j < n; j++){
            if (fabs(matrix[i][j]) > max_val ){
                max_val = fabs(matrix[i][j]);
                piv[0] = i;
                piv[1] = j;
            }
        }
    }
}


int sign(double num){
    if(num >= 0){ return 1; } else { return -1; }
}


double off(double ** matrix, int n){
    double sum = 0.0;
    int i,j;
    for (i=0 ; i < n; i++){
        for (j=0 ; j < n; j++){
            if (i != j){
                sum = sum + pow(matrix[i][j], 2);
            }
        }
    }
    return sum;
}


/* Assume that matrix argument is symetric */
void iter_jacobi(double ** matrix, double ** matrix_tag, double **eigenvectors, 
                        double **eigenvectors_tag, int n, int i, int j){
    int r;
    double t,c,s;
    double theta;
    
    theta = (matrix[j][j] - matrix[i][i]) / (2*matrix[i][j]);
    t = (sign(theta)) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    c = 1 / sqrt(pow(t, 2) + 1);
    s = t*c;

    /*printf("t= %f, c= %f, s= %f\n\n", t,c,s);*/

    /* Update Matrix Tag */
    for (r = 0; r < n; r++){
        if (r != i && r != j){
            matrix_tag[r][i] = c*matrix[r][i] - s*matrix[r][j];
            matrix_tag[r][j] = c*matrix[r][j] + s*matrix[r][i];
            matrix_tag[i][r] = matrix_tag[r][i];
            matrix_tag[j][r] = matrix_tag[r][j];
        }
    }
    matrix_tag[i][i] = pow(c, 2)*matrix[i][i] + pow(s, 2)*matrix[j][j] - 2*s*c*matrix[i][j];
    matrix_tag[j][j] = pow(s, 2)*matrix[i][i] + pow(c, 2)*matrix[j][j] + 2*s*c*matrix[i][j];
    matrix_tag[i][j] = 0;
    matrix_tag[j][i] = 0;

    /* Update Eigenvectors Matrix Tag */
    for (r = 0; r < n; r++) {
        eigenvectors_tag[r][i] = c * eigenvectors[r][i] - s * eigenvectors[r][j];
        eigenvectors_tag[r][j] = c * eigenvectors[r][j] + s * eigenvectors[r][i];
    }

    /* Update Eigenvectors Matrix */
    for (r = 0; r < n; r++) {
        eigenvectors[r][i] = eigenvectors_tag[r][i];
        eigenvectors[r][j] = eigenvectors_tag[r][j];
    }
    /*
    printf("eigen vectors:\n");
    print_mat(eigenvectors,n,n);
    printf("\ntag: \n");
    print_mat(eigenvectors_tag, n,n);*/
}


void update_matrix(double ** matrix, double ** matrix_tag, int n,
                                                        int i, int j){
    int r;
    for (r = 0; r < n; r++){
        matrix[r][i] = matrix_tag[r][i];
        matrix[r][j] = matrix_tag[r][j];
        matrix[i][r] = matrix_tag[r][i];
        matrix[j][r] = matrix_tag[r][j];
    }
}


int check_convergence(double ** matrix, double ** matrix_tag, int n){
    double epsilon = 1.0*pow(10,-5);
    if (off(matrix, n) - off(matrix_tag, n) <= epsilon){
        return 1;
    }
    return 0;
}


/* Assume that matrix argument is symetric */
void jacobi(double ** matrix, int n, double ** eigenvectors,
                                                double * eigenvalues){
    double ** matrix_tag;
    double ** eigenvectors_tag;
    int *piv;
    int i,j,r,k;
    int count_iter = 1;

    /* Find Pivot */
    piv = calloc(2, sizeof(int));
    assertion_check(piv == NULL);

    find_largest_element_pivot(matrix, n, piv);
    i = piv[0];
    j = piv[1];

    /* Initialize Matrix Tag */
    matrix_tag = make_mat(n, n);
    for (r = 0; r < n ; r++){
        for (k = 0; k < n ; k++){
            matrix_tag[r][k] = matrix[r][k];
        }
    }

    /* Initialize Eigenvectors Matrix Tag */
    eigenvectors_tag = make_mat_identity(n, n);

    /* Update Matrix Tag & Eigenvectors Matrix Tag & Eigenvectors Matrix */
    iter_jacobi(matrix, matrix_tag, eigenvectors, eigenvectors_tag, n, i, j);

    
    while (check_convergence(matrix, matrix_tag, n) != 1 
                            && count_iter <= MAX_ITER_JACOBI){
       /* printf("iter num: %d\n", count_iter);
        print_mat(eigenvectors, n, n);
        printf("\n\n");*/

        /* Update Matrix */
        update_matrix(matrix, matrix_tag, n, i, j);

        /* Update Pivot */
        find_largest_element_pivot(matrix, n, piv);
        i = piv[0];
        j = piv[1];

        /* Update Matrix Tag & Eigenvectors Matrix Tag & Eigenvectors Matrix */
        iter_jacobi(matrix, matrix_tag, eigenvectors, eigenvectors_tag, n,i,j);

        count_iter++;
    }

    /* Matrix Tag is Diagonal Matrix */
    for (r = 0; r < n; r++) {
        eigenvalues[r] = matrix_tag[r][r];
    }

    free(piv);
    free(matrix_tag);
    free(eigenvectors_tag);
}


void sort_descending(double * array, int n){
    int i, j;
    double temp_value;
    for (i = 0; i < n; ++i)
    {
        for (j = i+1; j < n; ++j)
        {
            if (array[i] < array[j])
            {
                temp_value = array[i];
                array[i] = array[j];
                array[j] = temp_value;
            }
        }
    }
}

/* assuming that eigenvalues are in descending order */
int largest_k_eigenvectors(double * eigenvalues, int n){
    int i, max_delta, k, N;
    double delta;

    max_delta = -1;
    k = 0;
    N = n/2;

    for(i=0; i<N; i++){
        delta = abs(eigenvalues[i]- eigenvalues[i+1]);
        if( delta > max_delta ){
            max_delta = delta;
            k = i + 1;
        }
    }
    return k;
}


void make_U(double ** eigenvectors, double * eigenvalues, 
                    double * eigenvalues_sorted, double ** U, int N, int K){
    int i,j,r;
    for (i=0 ; i < K ; i++){
        for (r = 0 ; r < N ; r++){
            if (eigenvalues_sorted[i] == eigenvalues[r]){
                for (j=0 ; j<N ; j++){
                    U[j][i] = eigenvectors[j][r];
                }
                eigenvalues[r] = -1;
                break;
            }
        }
    }
}


void make_T(double ** U, double ** T, int N, int K){
    int i,j;
    double * zero_vector = make_vector(N);
    double norm;
    for (i=0 ; i<N ; i++){
        norm = distance(U[i], zero_vector, K);
        for (j=0 ; j<K ; j++){
            T[i][j] = U[i][j]/norm;
        }
    }
    free(zero_vector);
}


/* KMEANS FUNCTIONS FOR PYTHON FILE */

/* allocate space for this array */
int * make_cluster_count(int K){
    int *count, i;
    count = calloc(K, sizeof(int));
    assertion_check(count == NULL);

    for(i = 0 ; i < K ; i++){
        count[i] = 0;
    }
    return count;
}


/* return the index of the closest cluster to this vector */
int get_closest_cluster(double *vector, double **centroids_array,
                                            int vector_size, int K){
    int closest_cluster = -1, j;
    double min_dist = DBL_MAX, dist;
    for(j=0 ; j < K ; j++){
        dist = distance(vector, centroids_array[j], vector_size);
        if(min_dist > dist){
            min_dist = dist;
            closest_cluster = j;
        }
    }
    return closest_cluster;
}


/* sums the vectors that in a certain cluster to this point of time */
void add_vector_to_cluster(double *vector, double *cluster_sum,
                                                    int vector_size){
    int j;
    for (j=0 ; j < vector_size ; j++){
        cluster_sum[j] += vector[j];
    }
}


/* use the information from cluster_sum array and cluster_count array to calculate the new centroids for each cluster */
int update_medians(double **cluster_sum, double **centroids_array,
                        int *cluster_count, int vector_size, int K){
    double epsilon = 0, norma, median;
    int epsilon_bool = 0;
    int i,j;
    for(i=0 ; i < K ; i++){
            norma = 0.0;
            for (j=0 ; j < vector_size ; j++){
                /*DRAFT: if(cluster_count[i] == 0){printf("An Error Has Occurred"); return 0;}*/
                assertion_check(cluster_count[i] == 0);
                median = cluster_sum[i][j]/cluster_count[i];
                norma += (centroids_array[i][j] - median)
                        *(centroids_array[i][j] - median);
                centroids_array[i][j] = median;

                /* resets the sum and count arrays for the next iteration */
                cluster_sum[i][j] = 0.0;
            }
            cluster_count[i] = 0;

            /* if there at least one vector that not close enough to median by epsilon */
            if(sqrt(norma) >= epsilon){
                epsilon_bool = 1;
            }
        }
    return epsilon_bool;
}


int k_means_c(double **vectors_array, double **centroids_array, int dimension,
                                    int num_of_vectors, int k, int max_iter){
    int iteration, i, closest_cluster;
    int *cluster_count;
    double **cluster_sum;

    cluster_sum = make_mat(k, dimension);
    cluster_count = make_cluster_count(k);

    if(cluster_sum == 0 || cluster_count == 0){
        return 0;
    }

    iteration = 0;
    while (iteration < max_iter){
        for (i=0 ; i < num_of_vectors ; i++){
            closest_cluster = get_closest_cluster(vectors_array[i], centroids_array, dimension, k);
            add_vector_to_cluster(vectors_array[i], cluster_sum[closest_cluster], dimension);
            cluster_count[closest_cluster]++;
        }

        /* also checks if we got under epsilon */
        if (!update_medians(cluster_sum, centroids_array, cluster_count, dimension, k)){
            break;
        }
        iteration++;
    }

    for(i=0 ; i < k ; i++){
        free(cluster_sum[i]);
    }
    free(cluster_sum);
    free(cluster_count);

    return 1;
}