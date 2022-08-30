#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "spkmeans.h"


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


static void make_jacobi_final_array(double ** eigenvectors, double * eigenvalues, double ** final_array, int n){
    int i,j;
    final_array = make_mat(n+1, n);

    for (i = 0 ; i < n ; i++){
        final_array[0][i] = eigenvalues[i];
    }
    for(i=0; i < n ; i++){
       for (j=0 ; j < n ; j++){
            final_array[i+1][j] = eigenvectors[i][j];
       }
    }
}


static PyObject* goal_fit(PyObject* self, PyObject* args){
    int dimension, n, i, goal_num, K, rows, culs;
    int rows_addition = 0;
    double **vectors_array, **final_array;
    double *eigenvalues = NULL;
    double *eigenvalues_sorted = NULL;
    double **eigenvectors = NULL;
    double **U = NULL;
    double **T = NULL;
    Graph * graph = malloc(sizeof(Graph*));
    Goal goal;

    PyObject * py_vectors_array;
    PyObject * py_final_array;

    /* get arguments to variables */
    if (!PyArg_ParseTuple(args, "iiiiO", &dimension, &n, &goal_num, &K, &py_vectors_array)){
        return NULL;
    }

    assertion_check(py_vectors_array == NULL);

    /* making c arrays with python input */
    vectors_array = convert_list_to_c(py_vectors_array, n, dimension);
    assertion_check(vectors_array == NULL);

    graph->dim = dimension;
    graph->n = n;
    graph->vertices = vectors_array;
    graph->ddg_mat = NULL;
    graph->lnorm_mat = NULL;
    graph->wam_mat = NULL;

    rows = n;
    culs = n;

    if (goal != e_jacobi){
        wam(graph);
        if (goal == e_wam){
            final_array = graph->wam_mat;
        }

        ddg(graph);
        if (goal == e_ddg){
            final_array = graph->ddg_mat;
        }

        lnorm(graph);
        if(goal == e_lnorm){
            final_array = graph->lnorm_mat;
        }

        if(goal == e_spk){
            eigenvectors = make_mat(n ,n);
            eigenvalues = make_vector(n);
            jacobi(graph->lnorm_mat, eigenvectors, eigenvalues, n);

            copy_mat(eigenvalues, eigenvalues_sorted, 1, n);
            sort_descending(eigenvalues_sorted, n);
            
            if (K == 0){
                K = largest_k_eigenvectors(eigenvalues_sorted, n);
            }

            U = make_mat(n, K);
            T = make_mat(n, K);
            make_U(eigenvectors, eigenvalues, eigenvalues_sorted, U, n, K);
            make_T(U, T, n, K);

            final_array = T;
            culs = K;
        }
    }

    else{
        eigenvectors = make_mat(n ,n);
        eigenvalues = make_vector(n);
        jacobi(graph->vertices, eigenvectors, eigenvalues, n);

        make_jacobi_final_array(eigenvectors, eigenvalues, final_array, n);
        rows = n+1;
    }

    py_final_array = convert_array_to_py(final_array, rows, culs);

    free_mat(vectors_array, n);
    free_mat(final_array, rows); /* what happens if we free the same thing twice? */
    free_mat(eigenvectors, n);
    free_mat(graph->vertices, n);
    free_mat(graph->wam_mat, n);
    free_mat(graph->ddg_mat, n);
    free_mat(graph->lnorm_mat, n);
    free_mat(U, n);
    free_mat(T, n);
    free(graph);
    free(eigenvalues);
    free(eigenvalues_sorted);
    
    return py_final_array;
}


static PyObject* fit_spk(PyObject* self, PyObject* args){
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
    {"goal_fit",
    (PyCFunction) goal_fit,
    METH_VARARGS,
    PyDoc_STR("if goal != spk")},
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