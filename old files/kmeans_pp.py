import sys
import numpy as np
import pandas as pd
import mykmeanssp as km


def error(b):
    sys.exit("Invalid Input!") if b else sys.exit("An Error Has Occurred")


def get_args():
    # checks if the arguments are valid
    n = len(sys.argv) - 5
    if n != 0 and n != 1:
        error(True)
        return
    try:
        num_of_clusters = int(sys.argv[1])
        max_iter = int(sys.argv[2]) if n == 1 else 300
        eps = float(sys.argv[2 + n])
        file_name_1 = str(sys.argv[3 + n])
        file_name_2 = str(sys.argv[4 + n])

        if num_of_clusters <= 0 or max_iter < 0 or eps < 0:
            error(True)

        return num_of_clusters, max_iter, eps, file_name_1, file_name_2

    except:
        error(True)


def create_inner_join(file1, file2):
    try:
        df_file1 = pd.read_csv(file1, header=None)
        df_file2 = pd.read_csv(file2, header=None)
        df_combine = pd.merge(df_file1, df_file2, on=0)
        df_combine = df_combine.sort_values(by=0)

        return df_combine.to_numpy()
    except:
        error(True)


def get_min_dist(vector, centroids_list, i):
    min_dist = float('inf')
    for j in range(i):
        temp = vector[1:] - centroids_list[j][1:]
        dist = np.dot(temp.T, temp)
        if min_dist > dist:
            min_dist = dist
    return min_dist


def k_means_plus_plus(vectors, num_of_clusters):
    # d - dimension, n - number of vectors
    probability = np.zeros(len(vectors))
    n = len(vectors)
    d = len(vectors[0])
    np.random.seed(0)
    cluster_centroids = np.zeros((num_of_clusters, d), float)
    index = np.random.choice(n)
    cluster_centroids[0] = vectors[index]
    for i in range(1, num_of_clusters):
        d_sum = 0
        for l in range(n):
            temp = get_min_dist(vectors[l], cluster_centroids, i)
            d_sum += temp
            probability[l] = temp
        probability = probability / d_sum
        index = np.random.choice(n, p=probability)
        cluster_centroids[i] = vectors[index]

    return cluster_centroids


def printSolution(final_centroids, observ_array):
    temp_txt = ""
    for i in range(len(observ_array)):
        temp_txt += str(observ_array[i][0])+","
    print(temp_txt[:-1])
    for i in range(len(final_centroids)):
        temp_txt = ""
        for j in range(len(final_centroids[0])):
            temp_txt += format(final_centroids[i][j], ".4f")+","
        print(temp_txt[:-1])


# write the numpy format
float_formatter = "{:.4f}".format
np.set_printoptions(formatter={'float_kind': float_formatter})

# checking arguments and putting them into variables
num_of_clusters, max_iter, eps, file_name_1, file_name_2 = get_args()
vectors = create_inner_join(file_name_1, file_name_2)

if num_of_clusters > len(vectors):
    error(True)

initial_centroids = k_means_plus_plus(vectors, num_of_clusters)
observ_index = initial_centroids[:, 0:1].astype('int32')
initial_centroids = initial_centroids[:, 1:].tolist()
vectors = vectors[:,1:].tolist()

final_centroids = km.fit(len(vectors[0]), len(vectors), num_of_clusters, max_iter, eps, vectors, initial_centroids)

printSolution(final_centroids, observ_index)