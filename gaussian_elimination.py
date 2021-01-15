# Ax = b
# A is the coefficient matrix.
# x is the unknown vector.
# b is the load vector.


def initialize():
    global answers, b_new, files, load_vector, matrices, a_coef
    answers = []                # answers to equations.(x vector)
    b_new = []                  # compare this to the original load vector to calculate numerical error.
    files = [open("A_1.txt", "r").read().split('\n'),
             open("A_2.txt", "r").read().split('\n'),
             open("A_3.txt", "r").read().split('\n')]

    # b vector
    load_vector = [list(map(float, open("b_1.txt", "r").read().split('\n'))),
                   list(map(float, open("b_2.txt", "r").read().split('\n'))),
                   list(map(float, open("b_3.txt", "r").read().split('\n')))]

    matrices = []   # coefficient matrices(A vectors) will be stored in this 3D array. each matrix will be augmented
                    # with corresponding load vector.
    a_coef = []     # coefficient matrices(A vectors). this array will not be augmented with load vector.


def read_augment_inputs():
    global answers, b_new, files, load_vector, matrices, a_coef
    for i in range(len(files)):
        matrices.append([])
        a_coef.append([])
        temp = []
        temp2 = []
        for j in range(len(files[i])):
            files[i][j] = [x for x in files[i][j].split(' ') if x != '']
            temp.append(list(map(float, files[i][j])))
            temp2.append(temp[j][:])
            temp[j].append(load_vector[i][j])                   # augment the matrix with load vector.
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


def scaled_partial_pivoting(matrix, col):
    scaling_vector = [abs(max(x, key=abs)) for x in matrix]
    pivot = abs(matrix[col][col]/scaling_vector[col])
    index = col
    for i in range(col, len(matrix)):
        if abs(matrix[i][col]/scaling_vector[i]) > pivot:
            index = i
            pivot = matrix[i][col]/scaling_vector[i]
    temp = matrix[col]
    matrix[col] = matrix[index]
    matrix[index] = temp


def gaussian_elimination(matrix, pivoting_function):
    assert len(matrix) == (len(matrix[0]) - 1)                  # check the size of augmented matrix.
    try:
        matrix_size = len(matrix)                   # matrix_size == number of rows.
        # forward elimination
        for col in range(matrix_size-1):
            pivoting_function(matrix, col)
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


def gauss_jordan(matrix, pivoting_function):
    assert len(matrix) == (len(matrix[0]) - 1)                  # check the size of augmented matrix.
    try:
        matrix_size = len(matrix)                   # matrix_size == number of rows.
        # forward elimination
        for col in range(matrix_size-1):
            pivoting_function(matrix, col)
            for row in range(col+1, matrix_size):
                multiple = -1 * (matrix[row][col] / matrix[col][col])
                matrix[row][col] = 0                    # explicitly assign zero to this element to avoid numerical errors.
                for k in range(col+1, matrix_size+1):                   # increase the range by 1 to include the load.
                    matrix[row][k] += multiple * matrix[col][k]
        # backward elimination
        for col in range(matrix_size-1, 0, -1):
            for row in range(col-1, -1, -1):
                multiple = -1 * (matrix[row][col] / matrix[col][col])
                for k in range(matrix_size, col-1, -1):
                    matrix[row][k] += multiple * matrix[col][k]
                matrix[row][col] = 0                    # explicitly assign zero to this element to avoid numerical errors.
        answer = [[0] for x in range(matrix_size)]
        for i in range(matrix_size):
            answer[i] = [matrix[i][matrix_size] / matrix[i][i]]
        return answer
    except Exception as e:
        if e.__class__ == ZeroDivisionError:
            print("Inconsistent equation")
        else:
            print("something else went wrong. ", e.__str__())


if __name__ == "__main__":
    global answers, b_new, files, load_vector, matrices, a_coef
    pivot_functions = [partial_pivoting, scaled_partial_pivoting]
    print_array = ["Gaussian elimination with partial pivoting :", "Gaussian elimination with scaled partial pivoting"]
    # gaussian elimination with different pivoting functions
    for k in range(2):
        print(print_array[k])
        initialize()
        read_augment_inputs()
        for i in range(3):
            answers.append(gaussian_elimination(matrices[i], pivot_functions[k]))
            b_new.append(matrix_multiplication(a_coef[i], answers[i]))
            error_norm = 0.0
            for j in range(len(load_vector[i])):
                error_norm += abs(load_vector[i][j] - b_new[i][j][0])
            for j in range(len(answers[i])):
                print("x"+str(j), " : ", answers[i][j])
            print("l1-norm of error vector : ", error_norm,"\n")
        print(100 * "/")
        # gauss-jordan method
        print("Gauss-jordan method")
        initialize()
        read_augment_inputs()
        for i in range(3):
            answers.append(gauss_jordan(matrices[i], partial_pivoting))
            b_new.append(matrix_multiplication(a_coef[i], answers[i]))
            error_norm = 0.0
            for j in range(len(load_vector[i])):
                error_norm += abs(load_vector[i][j] - b_new[i][j][0])
            for j in range(len(answers[i])):
                print("x" + str(j), " : ", answers[i][j])
            print("l1-norm of error vector : ", error_norm, "\n")
        print(100 * "/")


