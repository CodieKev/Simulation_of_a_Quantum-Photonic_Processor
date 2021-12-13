#!/usr/bin/env python
# coding: utf-8

# In[72]:


from math import *
import numpy as np
from itertools import product,permutations
from scipy.integrate import quad
from scipy import integrate
import matplotlib.pyplot as plt
R = 1/sqrt(2)
T = 1j/sqrt(2)
d = 1


# In[8]:


s = np.matrix([[R,T,0],[R*T,R**2,T],[T**2,R*T,R]])
u = np.transpose(s)


# In[44]:


def transform_multi(u):
    comb = list(product(range(len(u)), repeat=len(u)))  
    temp = []
    temp2 = []
    for i in range(len(comb)):
        temp3 = 1
        for j in range(len(comb[i])):
            temp3 = temp3*u[j,comb[i][j]]
        if temp3 !=0:
            temp.append(list(comb[i]))
            temp2.append(temp3)
    return(temp,temp2)
def id_state(temp1,pstate):
    temp2 = []
    ptemp = []
    l = len(temp1)
    while l >0:
        temp3 = temp1[0]
        temp4 = []
        temp5 =[]
        for j in range(len(temp1)):
            if sorted(temp1[j]) == sorted(temp3):
                temp4.append(temp1[j])
                temp5.append(j)
        temp2.append(temp4)
        temp1 = [x for x in temp1 if (x not in temp4)]
        ptemp.append([pstate[x] for x in temp5])
        pstate = [pstate[x] for x in range(len(pstate)) if (x not in temp5)]
        l = len(temp1)
    return(temp2,ptemp)
def time_list(s):
    tim = [list(range(len(s[0])))]*len(s)
    timn = []
    for i in range(len(s)):
        temp = s[i]
        temp1 = s[0]
        temp2 = tim[0]
        temp3 = []
        for j in range(len(temp)): 
            ind = temp1.index(temp[j])
            temp3.append(temp2[ind])
            temp1 = [temp1[x] for x in range(len(temp1)) if (x is not ind)]
            temp2 = [temp2[x] for x in range(len(temp2)) if (x is not ind)]
        timn.append(temp3)
    return(timn)
def replication(se):
    if len(se) - len(set(se)) !=0:
        temp = list(permutations(list(range(len(se)))))
        temp1 = []
        for i in range(len(temp)):
            temp2 = []
            for j in range(len(temp[i])):
                temp2.append(se[temp[i][j]])
            if temp2 == se:
                temp1.append(temp[i])
    else:
        temp1 = [se]
    return(temp1)
def parameters(rep,t,p):
    temp = []
    ptemp = []
    for i in range(len(rep)):
        temp2 = rep[i]
        for j in range(len(t)):
            temp3 = t[j]
            temp4 = [0]*len(temp3)
            for k in range(len(temp2)):
                if k == temp2[k]:
                    temp4[temp3.index(k)] = k
                else:
                    temp4[temp3.index(k)] = temp2[k]
            temp.append(temp4)
            ptemp.append(p[j]/len(rep))
    return(temp,ptemp)
def state_parameters(state,sp,specific_state):
    s1,p1 = id_state(state,sp)
    ind = 0
    for i in range(len(s1)):
        if sorted(specific_state) == sorted(s1[i][0]):
            ind = i
    s = s1[ind]
    p = p1[ind]
    t = time_list(s)
    ind_state = s[0]
    rep = replication(ind_state)
    par1,par2 = parameters(rep,t,p)
    return(par1,par2)


# In[45]:


state,sp = transform_multi(u)


# In[46]:


time_state,V = state_parameters(state,sp,[0,1,2])


# # Beta

# In[94]:


def epsl(t,a):
    result = (2*d**2/pi)**(1/4)*exp(-d**2*(a-t)**2)
    return(result)
def beta(x,y,z,t1,t2,t3):
    result = epsl(x,t1)*epsl(y,t2)*epsl(z,t3)
    return(result)
def gama(x,y,z):
    add = 0
    for i in range(len(V)):
        add += V[i]*beta(x,y,z,Stamp[3*i+0],Stamp[3*i+1],Stamp[3*i+2])
    result = (np.conj(add)*(add)).real
    return(result)
def time_stamp(time_state,time_mode):
    temp = []
    for i in range(len(time_state)):
        for j in range(len(time_state[i])):
            temp.append(time_mode[time_state[i].index(j)])
    return(temp)


# In[99]:


Tim = np.linspace(-4,4,51)


# In[101]:


p123 = [0]*51
time_state,V = state_parameters(state,sp,[0,1,2])
for i in range(51):
    print(i)
    t = Tim[i]
    time_mode = [0,t,8]
    Stamp = time_stamp(time_state,time_mode)
    prob = integrate.tplquad(gama, -8, 12, lambda x: -8, lambda x: 12,lambda x, y: -8, lambda x, y: 12)[0]
    p123[i] = prob


# In[102]:


p133 = [0]*51
time_state,V = state_parameters(state,sp,[0,2,2])
for i in range(51):
    print(i)
    t = Tim[i]
    time_mode = [0,t,8]
    Stamp = time_stamp(time_state,time_mode)
    prob = integrate.tplquad(gama, -8, 12, lambda x: -8, lambda x: 12,lambda x, y: -8, lambda x, y: 12)[0]
    p133[i] = prob


# In[103]:


p122 = [0]*51
time_state,V = state_parameters(state,sp,[0,1,1])
for i in range(51):
    print(i)
    t = Tim[i]
    time_mode = [0,t,8]
    Stamp = time_stamp(time_state,time_mode)
    prob = integrate.tplquad(gama, -8, 12, lambda x: -8, lambda x: 12,lambda x, y: -8, lambda x, y: 12)[0]
    p122[i] = prob


# In[104]:


Hom = np.array(p133)+np.array(p123)+np.array(p122)


# In[105]:


plt.subplot(221)
plt.plot(Tim, p122)
plt.subplot(222)
plt.plot(Tim, p123)
plt.subplot(223)
plt.plot(Tim, p133)
plt.subplot(224)
plt.plot(Tim, Hom)
plt.show()


# In[ ]:




