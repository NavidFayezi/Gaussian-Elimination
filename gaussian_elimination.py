files = [open("A_1.txt", "r").read().split('\n'),
         open("A_2.txt", "r").read().split('\n'),
         open("A_3.txt", "r").read().split('\n')]

load_vector = [list(map(float, open("b_1.txt", "r").read().split('\n'))),
               list(map(float, open("b_2.txt", "r").read().split('\n'))),
               list(map(float, open("b_3.txt", "r").read().split('\n')))]

matrices = []                   # coefficient matrices will be stored in this 3D array.
a_coef = []                     # coefficient matrices. this array will not be augmented with load vector.
for i in range(len(files)):
    matrices.append([])
    a_coef.append([])
    temp = []
    temp2 = []
    for j in range(len(files[i])):
        files[i][j] = [x for x in files[i][j].split(' ') if x != '']
        temp.append(list(map(float, files[i][j])))
        temp2.append(temp[j][:])
        temp[j].append(load_vector[i][j])                   # augment the matrix with load.
    matrices[i] = temp[:]
    a_coef[i] = temp2[:]


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


def gaussian_elimination_partial_pivoting(matrix):
    assert len(matrix) == (len(matrix[0]) - 1)                  # check the size of augmented matrix.
    try:
        matrix_size = len(matrix)                   # matrix_size == number of rows.
        # forward elimination
        for col in range(matrix_size-1):
            partial_pivoting(matrix, col)
            for row in range(col+1, matrix_size):
                multiple = -1 * (matrix[row][col] / matrix[col][col])
                matrix[row][col] = 0                    # explicitly assign zero to this element to avoid numerical errors.
                for k in range(col+1, matrix_size+1):                   # increase the range by 1 to include the load.
                    matrix[row][k] += multiple * matrix[col][k]
        # backward substitution
        answer = [[0] for x in range(matrix_size)]
        answer[matrix_size-1] = [(matrix[matrix_size-1][matrix_size] / matrix[matrix_size-1][matrix_size-1])]
        for i in range(matrix_size-2, -1, -1):
            sigma = 0.0
            for j in range(i+1, matrix_size):
                sigma += matrix[i][j] * answer[j][0]
            answer[i] = [((matrix[i][-1] - sigma) / matrix[i][i])]
        return answer
    except Exception as e:
        if e.__class__ == ZeroDivisionError:
            print("Inconsistent equation")
        else:
            print("something else went wrong. ", e.__str__())


def matrix_multiplication(a, b):                    # a * b
    assert len(a[0]) == len(b)
    result = []
    row_no = len(a)
    col_no = len(b[0])
    inner_loop = len(a[0])
    for row in range(row_no):
        result.append([])
        for col in range(col_no):
            temp = 0.0
            for k in range(inner_loop):
                temp += a[row][k] * b[k][col]
            result[row].append(temp)
    return result


a1_answer = gaussian_elimination_partial_pivoting(matrices[0])
b1_new = matrix_multiplication(a_coef[0], a1_answer)
print(b_new)