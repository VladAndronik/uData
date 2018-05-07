import copy
import numpy as np

class GaussSolver(object):
    
    __matrix = None
    __b = None
    __solution = []
    __err = None
    
    
    def __init__(self, matrix, b):
        self.__matrix = matrix
        self.__b = b
    
    def _concatenate(self, A, b):
        i = 0
        for el in A:
            el.append(b[i])
            i += 1
        return A
    
    
    def toString(self, A):
        for a in A:
            print('\t','|', a, '|', '\n')
                
    
    def linEqSolver(self):
        
        """
            Solve the linear equation Ax = b
            Parameters
            ------
            __matrix : list of lists 
                       Coefficients matrix
            __b : list 
                  Right-hand vector
            
            Returns
            ------
            X : list of equation solutions
            error : integer
                    code of the error
                    -1 : error while executing the loop
                    0  : everything all right
        """
        
        error = 0
        M = copy.deepcopy(self.__matrix)
        b = copy.deepcopy(self.__b)
        A = self._concatenate(M, b)
        
        n = len(A)
        m = len(A[0])
        
    
        for i in range(0, n):
            
            # Find max value in the current column
            maxEl = abs(A[i][i])
            maxRow = i
            for k in range(i+1, n):
                if abs(A[k][i]) > maxEl:
                    maxEl = abs(A[k][i])
                    maxRow = k

            # Swap max row with current row
            for k in range(i, n+1):
                temp = A[maxRow][k]
                A[maxRow][k] = A[i][k]
                A[i][k] = temp
            
            # Zero out the rows below the current column
            for k in range(i+1, n):
                if A[i][i] == 0:
                    error = -i
                    break
                c = -A[k][i]/A[i][i]
                for j in range(i, n+1):
                    if i == j:
                        A[k][j] = 0
                    else:
                        A[k][j] += c * A[i][j]
            print('Matrix of coefficients after %s step of forward pass: '%str(i+1))
            print('\t')
            self.toString(A)

        # Find the roots for upper-triangular matrix
        X = [0 for i in range(n)]
        for i in range(n-1, -1, -1):
            X[i] = A[i][n]/A[i][i]
            for k in range(i-1, -1, -1):
                A[k][n] -= A[k][i] * X[i]
        self.__solution = X
        self.__err = error
            
    def verification(self, A, b):
        X = np.linalg.solve(A, b)
        return X
    
    def _r(self, A, x, b):
        r = np.matmul(np.array(A), np.array(x)) - np.array(b)
        return r
    
    def toStringAll(self):
        
        print('Matrix of coefficients:\n')
        self.toString(self.__matrix)
        print('Right-hand vector: \n')
        self.toString(self.__b)
        
        print('Vector of solutions with Gauss class: \n')
        self.toString(self.__solution)
        print('Code of error: %i'%self.__err)
        
        print('\nVerification with NumPy: \n')
        self.toString(self.verification(np.array(self.__matrix), np.array(self.__b)))
        
        print('\nCheck Ax - b = r \n')
        self.toString(self._r(self.__matrix, self.__solution, self.__b))
    
    