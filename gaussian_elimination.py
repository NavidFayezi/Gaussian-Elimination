files = []
files.append(open("A_1.txt", "r").read().split('\n'))
files.append(open("A_2.txt", "r").read().split('\n'))
files.append(open("A_3.txt", "r").read().split('\n'))

matrices = []   # all three matrices will be stored in this 3D array
for i in range(len(files)):
    matrices.append([])
    temp = []
    for j in range(len(files[i])):
        files[i][j] = [x for x in files[i][j].split(' ') if x != '']
        temp.append(list(map(float, files[i][j])))
    matrices[i] = temp[:]


def partial_pivoting(matrix, col):
    pivot = matrix[col][col]
    index = col
    for i in range(col, len(matrix)):
        if matrix[i][col] > pivot:
            index = i
            pivot = matrix[i][col]
    temp = matrix[col]
    matrix[col] = matrix[index]
    matrix[index] = temp


