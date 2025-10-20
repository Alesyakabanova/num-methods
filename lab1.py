import numpy as np
import time

def get_coeff(a, b, c):
    n = len(b)
    diag = np.zeros(n)
    alpha = np.zeros(n-1)
    
    diag[0] = b[0]
    for i in range(n-1):
        alpha[i] = c[i] / diag[i]       
        diag[i+1] = b[i+1] - a[i] * alpha[i]  
    
    detA = np.prod(diag)
    return diag, alpha, detA


def find_one_column(a, alpha, diag, f):
   
    n = len(f)
    y = np.zeros(n)
    y[0] = f[0]
    for i in range(1, n):
        y[i] = f[i] - a[i-1] * y[i-1] / diag[i-1]  
    
    x = np.zeros(n)
    x[-1] = y[-1] / diag[-1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - alpha[i] * x[i+1]) / diag[i] 
    
    return x


def find_inverse(a, b, c):
    
    n = len(b)
    diag, alpha, detA = get_coeff(a, b, c)
    invA = np.zeros((n, n))
    
    for i in range(n):
        f = np.zeros(n)
        f[i] = 1.0
        x = find_one_column(a, alpha, diag, f)
        invA[:, i] = x
    
    return invA, detA

c = np.array([1, 5, -2, 0, 3, 0, 3, 5, -3, -3, 3, 0, 4, 1, -1, -1]) 
d = np.array([5, 3, 6, -6, 6, -2, 3, 5, 6, -2, -6, 5, 2, 4, 1, -2, 4])  
e = np.array([-3, 2, -1, 3, -4, 0, -5, -1, 3, 2, -1, -4, -1, 5, -1, -2]) 

'''c = np.array([-4, 2, -2, 5, -5, -3, 0, 0, 5])
d = np.array([4, -6, 4, 3, -2, 1, -6, 2, 3, -4])
e = np.array([-3, 4, -4, -1, 4, -4, 0, 4, 2])'''


'''c = np.array([2, -2, -2])
d = np.array([-2, 3, 4, 3])
e = np.array([2, 5, 1])'''

invA, detA = find_inverse(c, d, e)

print("Определитель матрицы:")
print(detA)
print("Обратная матрица:")
print(invA)

def tridiag_example(n):
    a = np.random.uniform(-5, 5, n-1)
    b = np.random.uniform(-5, 5, n)
    c = np.random.uniform(-5, 5, n-1)
    return a, b, c

sizes = [50, 100, 200, 400, 800]
times = []

for n in sizes:
    a, b, c = tridiag_example(n)
    start = time.time()
    invA, detA = find_inverse(a, b, c)
    end = time.time()
    times.append(end - start)
    print(f"n={n:4d} | время={end - start:.5f} сек | det(A)={detA:.3e}")

print("\nРезультаты эксперимента:")
for n, t in zip(sizes, times):
    print(f"n={n:<5}  t={t:.5f} сек")



