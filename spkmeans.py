import sys
import numpy as np
import pandas as pd
import mykmeanssp as km


MAX_ITER = 300


def error(b):
    sys.exit("Invalid Input!") if b else sys.exit("An Error Has Occurred")


def get_args():
    # checks if the arguments are valid
    argc = len(sys.argv)
    if argc != 4:
        error(True)
        return
    try:
        k = int(sys.argv[1])
        goal_str = str(sys.argv[2])
        input_name = str(sys.argv[3])
        if goal_str != "spk" and goal_str != "wam" and goal_str != "ddg" and goal_str != "lnorm" and goal_str != "jacobi":
            error(True)
            return

        return k, goal_str, input_name

    except:
        error(True)


def printSolution_spk(final_centroids, observ_array):
    temp_txt = ""
    for i in range(len(observ_array)):
        temp_txt += str(observ_array[i][0])+","
    print(temp_txt[:-1])
    for i in range(len(final_centroids)):
        temp_txt = ""
        for j in range(len(final_centroids[0])):
            temp_txt += format(final_centroids[i][j], ".4f")+","
        print(temp_txt[:-1])


def printSolution(mat):
    temp_txt = ""
    for i in range(len(mat)):
        temp_txt = ""
        for j in range(len(mat[0])):
            temp_txt += format(mat[i][j], ".4f")+","
        print(temp_txt[:-1])


def get_min_dist(vector, centroids_list, i):
    min_dist = float('inf')
    for j in range(i):
        temp = vector[1:] - centroids_list[j][1:]
        dist = np.dot(temp.T, temp)
        if min_dist > dist:
            min_dist = dist
    return min_dist


def k_means_plus_plus(vectors, k):
    # d - dimension, n - number of vectors
    probability = np.zeros(len(vectors))
    n = len(vectors)
    d = len(vectors[0])
    np.random.seed(0)
    cluster_centroids = np.zeros((k, d), float)
    index = np.random.choice(n)
    cluster_centroids[0] = vectors[index]
    for i in range(1, k):
        d_sum = 0
        for l in range(n):
            temp = get_min_dist(vectors[l], cluster_centroids, i)
            d_sum += temp
            probability[l] = temp
        probability = probability / d_sum
        index = np.random.choice(n, p=probability)
        cluster_centroids[i] = vectors[index]

    return cluster_centroids


if __name__ == '__main__':
    # write the numpy format
    float_formatter = "{:.4f}".format
    np.set_printoptions(formatter={'float_kind': float_formatter})

    # checking arguments and putting them into variables
    k, goal_str, input_name = get_args()

    # getting the vectors from file
    try:
        input_vectors = pd.read_csv(input_name, header=None)
        input_vectors = pd.DataFrame(input_vectors).values.tolist()
        # input_vectors = input_vectors.to_numpy().tolist()
    except:
        error(True)

    dim = len(input_vectors[0])
    N = len(input_vectors)

    if k >= N or k < 0:
        error(True)

    final_solution_mat = km.goal_fit(dim, N, ord(goal_str[0]), k, input_vectors)

    k = len(final_solution_mat[0])

    # adding index cull
    for i in range(len(final_solution_mat)):
        final_solution_mat[i].insert(0, i)

    final_solution_mat = np.asarray(final_solution_mat)
    initial_centroids = k_means_plus_plus(final_solution_mat, k)
    observ_index = initial_centroids[:, 0:1].astype('int32')
    initial_centroids = initial_centroids[:, 1:].tolist()
    vectors = final_solution_mat[:, 1:].tolist()

    final_centroids = km.fit_kmeans(len(vectors[0]), len(vectors), k, MAX_ITER, vectors, initial_centroids)

    printSolution_spk(final_centroids, observ_index)
