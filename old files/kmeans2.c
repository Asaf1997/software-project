#define PY_SSIZE_T_CLEAN
#include <Python.h>
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

/* allocate space for this array */
static double ** make_cluster_sum(int K, int vector_size){
    double **cluster_sum;
    int i,j;
    cluster_sum = calloc(K+1, sizeof(double*));
    for (i = 0; i < K; i++){
        cluster_sum[i] = calloc(vector_size+1, sizeof(double));
        for (j=0 ; j < vector_size ; j++){
            cluster_sum[i][j] = 0.0;
        }
    }
    return cluster_sum;
}

/* allocate space for this array */
static int * make_cluster_count(int K){
    int *count, i;
    count = calloc(K, sizeof(int));
    for(i = 0 ; i < K ; i++){
        count[i] = 0;
    }
    return count;
}

/* return the index of the closest cluster to this vector */
static int get_closest_cluster(double *vector, double **centroids_array, int vector_size, int K){
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
static void add_vector_to_cluster(double *vector, double *cluster_sum, int vector_size){
    int j;
    for (j=0 ; j < vector_size ; j++){
        cluster_sum[j] += vector[j];
    }
}

/* use the information from cluster_sum array and cluster_count array to calculate the new centroids for each cluster */
static int update_medians(double **cluster_sum, double **centroids_array, int *cluster_count, int vector_size, int K){
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


static int k_means_c(double **vectors_array, double **centroids_array, int dimension, int num_of_vectors, int k, int max_iter, double epsilon){
    int iteration, i, closest_cluster;
    int *cluster_count;
    double **cluster_sum;

    cluster_sum = make_cluster_sum(k, dimension);
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


static double** convert_list_to_c(PyObject* pyList, int rows, int cull){
    double **array = calloc(rows, sizeof(double*));
    int index_i, index_j;
    PyObject * py_row;
    PyObject * py_num;


    if(!PyList_Check(pyList)){
        return NULL;
    }

    for (index_i = 0 ; index_i < rows ; index_i++){
        array[index_i] = calloc(cull, sizeof(double));
        py_row = PyList_GetItem(pyList, index_i);

        if (!PyList_Check(py_row)){
            return NULL;
        }

        for(index_j=0 ; index_j < cull ; index_j++){
            py_num = PyList_GetItem(py_row, index_j);

            if (!PyFloat_Check(py_num)){
                return NULL;
            }

            array[index_i][index_j] = PyFloat_AsDouble(py_num);
        }
    }

    return array;
}


static PyObject* convert_array_to_py(double** array, int rows, int cull){
    int index_i, index_j;

    PyObject * py_table = PyList_New(rows);
    PyObject * py_row;
    PyObject * py_num;

    for (index_i=0 ; index_i<rows ; index_i++){
        py_row = PyList_New(cull);
        for(index_j=0 ; index_j < cull ; index_j++){
            py_num = Py_BuildValue("d", array[index_i][index_j]);
            PyList_SetItem(py_row, index_j, py_num);
        }
        PyList_SetItem(py_table, index_i, py_row);
    }

    return py_table;
}


static PyObject* fit(PyObject* self, PyObject* args){
    int dimension, num_of_vectors, k, max_iter, i;
    double epsilon;
    double **vectors_array, **centroids_array;

    PyObject * py_vectors_array;
    PyObject * py_centroids_array;
    PyObject * py_final_centroids;

    /* get arguments to variables */
    if (!PyArg_ParseTuple(args, "iiiidOO", &dimension, &num_of_vectors, &k, &max_iter, &epsilon, &py_vectors_array, &py_centroids_array)){
        return NULL;
    }

    /* making c arrays with python input */
    vectors_array = convert_list_to_c(py_vectors_array, num_of_vectors, dimension);
    centroids_array = convert_list_to_c(py_centroids_array, k, dimension);

    if (vectors_array == NULL || centroids_array == NULL){
        return NULL;
    }


    /* initialize kmeans algorithem */
    if(!k_means_c(vectors_array, centroids_array, dimension, num_of_vectors, k, max_iter, epsilon)){
        return NULL;
    }

    py_final_centroids = convert_array_to_py(centroids_array, k, dimension);

    for(i=0 ; i < num_of_vectors ; i++){
        free(vectors_array[i]);
    }
    for(i=0 ; i < k ; i++){
        free(centroids_array[i]);
    }
    free(vectors_array);
    free(centroids_array);

    return py_final_centroids;
}


static PyMethodDef mykmeanssp_Methods[] = {
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("k means func")},
    {NULL,NULL,0,NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    mykmeanssp_Methods
};


PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}