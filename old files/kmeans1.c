#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

double distance(double vector1[], double vector2[], int length){
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

double ** make_centroids_array(double **array, int K, int vector_size){
    int i,j;
    double **centroids_array;
    centroids_array = calloc(K, sizeof(double*));
    for (i = 0; i < K; i++){
        centroids_array[i] = calloc(vector_size, sizeof(double));
        for (j=0 ; j < vector_size ; j++){
            centroids_array[i][j] = array[i][j];
        }
    }
    return centroids_array;
}

double ** make_cluster_sum(int K, int vector_size){
    double **cluster_sum;
    int i,j;
    cluster_sum = calloc(K, sizeof(double*));
    for (i = 0; i < K; i++){
        cluster_sum[i] = calloc(vector_size, sizeof(double));
        for (j=0 ; j < vector_size ; j++){
            cluster_sum[i][j] = 0.0;
        }
    }
    return cluster_sum;
}

int * make_cluster_count(int K){
    int *count, i;
    count = calloc(K, sizeof(int));
    for(i = 0 ; i < K ; i++){
        count[i] = 0;
    }
    return count;
}

int write_on_file(double **array, int max_i, int max_j , char file_name[]){
    FILE *file;
    int i,j;
    
    if((file = fopen(file_name, "w")) == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    for (i=0 ; i<max_i ; i++){
        for(j=0 ; j<max_j ; j++){
            if (j < max_j-1){
                fprintf(file, "%.4f,", array[i][j]); }
            else{
                fprintf(file, "%.4f\n", array[i][j]);
            }
        }
    }
    fclose(file);
    return 0;
}

int get_closest_cluster(double *vector, double **centroids_array, int vector_size, int K){
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

void add_vector_to_cluster(double *vector, double *cluster_sum, int vector_size){
    int j;
    for (j=0 ; j < vector_size ; j++){
        cluster_sum[j] += vector[j];
    }
}

int update_medians(double **cluster_sum, double **centroids_array, int *cluster_count, int vector_size, int K){
    double epsilon = 0.001, norma, median;
    int epsilon_bool = 0;
    int i,j;
    for(i=0 ; i < K ; i++){
            norma = 0.0;
            for (j=0 ; j < vector_size ; j++){
                if(cluster_count[i] == 0){
                    printf("An Error Has Occurred");
                    return 0;
                }
                median = cluster_sum[i][j]/cluster_count[i];
                norma += (centroids_array[i][j] - median)*(centroids_array[i][j] - median);
                centroids_array[i][j] = median;
                cluster_sum[i][j] = 0.0;
            }
            cluster_count[i] = 0;
            if(sqrt(norma) >= epsilon){
                epsilon_bool = 1;
            }
        }
    return epsilon_bool;
}

int main(int argc, char *argv[]){
    FILE *input_file;
    int num_of_vectors, vector_size, *cluster_count;
    int i, K, max_iter, iteration, closest_cluster;
    char *input_name, *output_name;
    double **vectors_array, **centroids_array, **cluster_sum;
    
    /* checking correctness of the input */
    if(argc < 4 || argc > 5 ){        
        printf("Invalid Input!");
        return 1;
    }
    if( !(K = atoi(argv[1])) ){
        printf("Invalid Input!");
        return 1;
    }
    if(argc == 4){
        max_iter = 200;
        input_name = argv[2];
        output_name = argv[3];
    }
    else{
        input_name = argv[3];
        output_name = argv[4];
        if(!(max_iter = atoi(argv[2])) ){
            printf("Invalid Input!");
            return 1;
        }
    }
    if((input_file = fopen(input_name, "r")) == NULL){
        printf("Invalid Input!");
        return 1;
    }

    vector_size = get_vector_size(input_file);
    num_of_vectors = get_num_of_vectors(input_file);
    vectors_array = make_vectors_array(input_file, vector_size, num_of_vectors);
    centroids_array = make_centroids_array(vectors_array, K, vector_size);
    cluster_sum = make_cluster_sum(K, vector_size);
    cluster_count = make_cluster_count(K);

    if(vectors_array == 0 || centroids_array == 0 || cluster_sum == 0 || cluster_count == 0){
        printf("An Error Has Occurred");
        return 1;
    }

    iteration = 0;
    while (iteration < max_iter){
        for (i=0 ; i < num_of_vectors ; i++){
            closest_cluster = get_closest_cluster(vectors_array[i], centroids_array, vector_size, K);
            add_vector_to_cluster(vectors_array[i], cluster_sum[closest_cluster], vector_size);
            cluster_count[closest_cluster]++;
        }
        if (!update_medians(cluster_sum, centroids_array, cluster_count, vector_size, K)){
            break;
        }
        iteration++;
    }

    if(write_on_file(centroids_array, K, vector_size, output_name)){
        return 1;
    }
    fclose(input_file);
    for(i=0 ; i < num_of_vectors ; i++){
        free(vectors_array[i]);
    }
    for(i=0 ; i < K ; i++){
        free(centroids_array[i]);
        free(cluster_sum[i]);
    }
    free(vectors_array);
    free(centroids_array);
    free(cluster_sum);
    free(cluster_count);
    return 0;
}