#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pandas import read_excel
from printarray import printarray


# In[2]:


KDF = read_excel('stiffness_matrix.xlsx')
MDF = read_excel('mass_matrix.xlsx')
CDF = read_excel('damping_matrix.xlsx')


# In[3]:


K = KDF.values
M = MDF.values
C = CDF.values


# In[4]:


print('Upper left sub-half-matrix of K')
for ii in range(4):
    for jj in range(4):print('{:.15f}  '.format(K[ii,jj]),end="")

    print('')

print('\nLower right sub-half-matrix of K')
for ii in range(4,8):
    for jj in range(4,8):print('{:.15f}  '.format(K[ii,jj]),end="")

    print('')

print('\nUpper right sub-half-matrix of K')
for ii in range(4):
    for jj in range(4,8):print('{:.15f}  '.format(K[ii,jj]),end="")

    print('')


# In[5]:


shouldbezero8times8 = abs(K - K.T)
print('shouldbezero8times8:')
for ii in range(8):
    printarray(shouldbezero8times8[ii,:],'.10f')


# In[6]:


print('\nUpper left sub-half-matrix of M')
for ii in range(4):
    for jj in range(4):print('{:.15f}  '.format(M[ii,jj]),end="")

    print('')

print('\nLower right sub-half-matrix of M')
for ii in range(4,8):
    for jj in range(4,8):print('{:.15f}  '.format(M[ii,jj]),end="")

    print('')


print('\nUpper right sub-half-matrix of M')
for ii in range(4):
    for jj in range(4,8):print('{:.15f}  '.format(M[ii,jj]),end="")

    print('')


# In[7]:


shouldbezero8times8 = abs(M - M.T)
print('shouldbezero8times8:')
for ii in range(8):
    printarray(shouldbezero8times8[ii,:],'.10f')


# In[8]:


print('\nUpper left sub-half-matrix of C')
for ii in range(4):
    for jj in range(4):print('{:.15f}  '.format(C[ii,jj]),end="")

    print('')

print('\nLower right sub-half-matrix of C')
for ii in range(4,8):
    for jj in range(4,8):print('{:.15f}  '.format(C[ii,jj]),end="")

    print('')


print('\nUpper right sub-half-matrix of C')
for ii in range(4):
    for jj in range(4,8):print('{:.15f}  '.format(C[ii,jj]),end="")

    print('')


# In[9]:


shouldbezero8times8 = abs(C - C.T)
print('shouldbezero8times8:')
for ii in range(8):
    printarray(shouldbezero8times8[ii,:],'.10f')

import numpy as np
from scipy.linalg import eigh

# --- 下面是添加的计算部分 ---

# 1. Solve the Generalized Eigenvalue Problem: (K - w^2*M)v = 0
# 使用 scipy.linalg.eigh，它专门用于求解对称矩阵 (Symmetric Matrices) 的广义特征值问题
# eigenvalues (w^2) 会自动按从小到大排序
evals, evecs = eigh(K, M)

# 2. Calculate Natural Frequencies (rad/s)
# evals 里面存的是 frequencies 的平方 (omega^2)，所以要开根号
natural_frequencies = np.sqrt(evals)

print("\n-------------------------------------------")
print("Results for Homework 5")
print("-------------------------------------------")

# 打印所有的 Natural Frequencies (作业第1部分)
print("Undamped Natural Frequencies (rad/s):")
for i, freq in enumerate(natural_frequencies):
    print(f"Mode {i+1}: {freq:.4f} rad/s")

# 3. Calculate the ratio for the lowest natural frequency (作业第2部分)
# eigh 返回的结果通常已经排序，第0个就是 lowest frequency
lowest_freq_idx = 0
lowest_mode_shape = evecs[:, lowest_freq_idx]

# 获取第一个位移 (index 0) 和最后一个位移 (index 7，因为是8自由度)
disp_first = lowest_mode_shape[0]
disp_last = lowest_mode_shape[7]

ratio = disp_first / disp_last

print(f"\nLowest Natural Frequency is: {natural_frequencies[lowest_freq_idx]:.4f} rad/s")
print(f"Displacement Ratio (1st / 8th) for this mode: {ratio:.6f}")