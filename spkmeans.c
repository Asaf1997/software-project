#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* returns the distance between two vectors */
static double distance(double vector1[], double vector2[], int length){
    double sum=0.0;
    int i;
    for(i=0; i < length; i++){
        sum += (vector1[i] - vector2[i])*(vector1[i] - vector2[i]);
    }
    return sum;
}

int get_vector_size(FILE * file){
    char line[1000];
    int counter = 1;
    int i;
    fseek(file, 0, SEEK_SET); 
    /* get to the start of the file */
    if(fgets(line, 1000, file)){
        for (i = 0 ; i < 1000 ; i++){
            if (line[i] == ','){
                counter++;
            }
        }
    }
    return counter;
}

int get_num_of_vectors(FILE * file){
    char line[1000];
    int counter = 0;
    fseek(file, 0, SEEK_SET);
    /* get to the start of the file */
    while (fgets(line, 1000, file)){
       counter++;
    }
    return counter;
}

double ** make_vectors_array(FILE * file, int vector_size, int num_of_vectors){
    int i=0,j=0;
    double **vectors_array;
    char trash;
    fseek(file, 0, SEEK_SET);  /* get to the start of the file */
    vectors_array = calloc(num_of_vectors, sizeof(double*));
    for(i=0 ; i < num_of_vectors ; i++){
        vectors_array[i] = calloc(vector_size, sizeof(double));
        for(j=0 ; j<vector_size ; j++){
            fscanf(file, "%lf", &vectors_array[i][j]);
            fscanf(file, "%c", &trash);
        }
    }
    return vectors_array;
}

static double ** make_mat(int n){
    double ** wam_mat;
    int i,j;
    wam_mat = calloc(n, sizeof(double*));
    
    if (wam_mat == NULL){
            printf("An Error Has Occurred");
            return NULL;
    }

    for (i = 0; i < n; i++){
        wam_mat[i] = calloc(n, sizeof(double));
        
        if (wam_mat[i] == NULL){
            printf("An Error Has Occurred");
            return NULL;
        }

        for (j=0 ; j < n ; j++){
            wam_mat[i][j] = 0.0;
        }
    }
    return wam_mat;
}

static void print_mat(double ** mat, int i, int j){
    int x,y;
    for (x = 0 ; x < i ; x++){
        for (y = 0 ; y < j ; y++){
            if (y == j-1){
                printf("%.4f", mat[i][j]);
            }
            else{
                printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }
}

static int wam(double ** vectors, int d, int n, char * goal){
    double ** wam_mat = make_mat(n);
    double weight;
    int i,j;

    if (wam_mat == NULL){
        return 1;
    }

    for (i = 0; i < n; i++){
        for (j=i+1 ; j < n ; j++){
            weight = exp(-sqrt(distance(vectors[i], vectors[j], d))/2);
            wam_mat[i][j] = weight;
            wam_mat[j][i] = weight;
        }
    }

    if (goal == "wam"){
        print_mat(wam_mat, n ,n);
        for(i=0 ; i < n ; i++){
            free(wam_mat[i]);
        }
        free(wam_mat);
        return 0;
    }
    else{
        return ddg(vectors, wam_mat, d, n, goal);
    }
}

static int ddg(double ** vectors, double ** wam_mat, int d, int n, char * goal){
    double ** weighted_adj_matrix = wam_mat;
    double ** diagonal_matrix = calloc(n, sizeof(double*));    
    double sum_weights = 0;
    int i;

    if (diagonal_matrix == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    for(int i=0 ; i < n ; i++){
        diagonal_matrix[i] = calloc(n, sizeof(double));
        
        if (diagonal_matrix[i] == NULL){ printf("An Error Has Occurred"); return 1;}

        for(int j=0 ; j < n ; j++){
            sum_weights += weighted_adj_matrix[i][j];
            if (i != j){
                diagonal_matrix[i][j] = 0.0;
            }
        }
        diagonal_matrix[i][i] = sum_weights;
    }
    
    if (goal == "ddg"){
        print_mat(diagonal_matrix, n, n);
        
        for(i=0 ; i < n ; i++){
            free(wam_mat[i]);
        }
        free(wam_mat);
        
        for(i=0 ; i < n ; i++){
            free(diagonal_matrix[i]);
        }
        free(diagonal_matrix);

        return 0;
    }

    else{
        return lnorm(vectors, wam_mat, diagonal_matrix, d, n, goal);
    }
}

static double ** sqrt_diagonal_matrix(double ** diagonal_matrix, int n){
    double ** sqrt_matrix = calloc(n, sizeof(double*));

    if (sqrt_matrix == NULL){ printf("An Error Has Occurred"); return NULL;}

    for(int i=0 ; i < n ; i++){
        sqrt_matrix[i] = calloc(n, sizeof(double));

        if (sqrt_matrix[i] == NULL){ printf("An Error Has Occurred"); return NULL;}

        for(int j=0 ; j < n ; j++){
            if (i == j) { sqrt_matrix[i][j] = 1 / sqrt(diagonal_matrix[i][i]); }
            else { sqrt_matrix[i][j] = 0.0; }
        }
    }
    return sqrt_matrix;
}

static double ** lnorm(double ** vectors, double ** wam_mat, double ** ddg_mat, int d, int n, char * goal){
    double ** sqrt_D = sqrt_diagonal_matrix(ddg_mat, n);
    double ** lnorm_mat = make_mat(n);
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

    if(goal == "lnorm"){
        print_mat(lnorm_mat, n ,n);

        for(i=0 ; i < n ; i++){
            free(wam_mat[i]);
        }
        free(wam_mat);
        
        for(i=0 ; i < n ; i++){
            free(ddg_mat[i]);
        }
        free(ddg_mat);

        for(i=0 ; i < n ; i++){
            free(ddg_mat[i]);
        }
        free(ddg_mat);

        for(i=0 ; i < n ; i++){
            free(lnorm_mat[i]);
        }
        free(lnorm_mat);

        for(i=0 ; i < n ; i++){
            free(sqrt_D[i]);
        }
        free(sqrt_D);

        return 0;
    }

    else{
        return jacobi(vectors, wam_mat, ddg_mat, sqrt_D, lnorm_mat , d, n, goal);
    }
}

static double ** jacobi(double ** vectors, double ** wam_mat, double ** ddg_mat, double ** sqrt_D, double ** lnorm_mat, int d, int n, char * goal){}


int main(int argc, char *argv[]){
    FILE *input_file;
    int num_of_vectors, vector_size, *cluster_count;
    int i, K, max_iter, iteration, closest_cluster;
    char *input_name, *goal;
    double **vectors_array, **centroids_array, **cluster_sum;

    /* checking correctness of the input */
    if(argc != 3 ){        
        printf("Invalid Input!");
        return 1;
    }

    input_name = argv[1];
    goal = argv[2];
 
    if((input_file = fopen(input_name, "r")) == NULL){
        printf("Invalid Input!");
        return 1;
    }

    vector_size = get_vector_size(input_file);
    num_of_vectors = get_num_of_vectors(input_file);
    vectors_array = make_vectors_array(input_file, vector_size, num_of_vectors);
}