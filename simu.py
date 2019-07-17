#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 00:35:22 2019

@author: viviantian
"""


import numpy as np
from random import choice
import pdb
import sys
import time
import matplotlib.pyplot as plt
import os
import copy
start = time.time()
time_arr = []
sys.setrecursionlimit(10000) #递归深度设置为10000 
# 聚合物（引发剂）的数量 n_Polymer, NPC（单体）的数量n_NPC, AcH的数量n_AcH，NCA的数量n_NCA
# NPC：polymer：ACH=20:1:10

n_k4 = 0 #记录发生k4反应的次数
n_k6 = 0 #记录发生k6反应的次数

Polymer_list = []  #聚合物列表
kk = 1
# 初始化参数
T = 0
n_Polymer_init = int(1000*kk)
n_NPC_init = int(20000*kk)
n_Polymer = n_Polymer_init
n_NPC = n_NPC_init
n_AcH = int(10000*kk)
n_NCA = 0
n_Ac = 0
n_NPCremain = 0 #酸化后剩余的NPC
#n_NH3 = 0
n_phenol = 0 # 苯酚的数目 = 反应掉的NPC数目
n_side = 0 # 侧基的数目
n_HA = n_NPC+ n_AcH; # 酸的总量
n_activechain = int(1000*kk) #活性链的数目初始和引发剂的数目相同
kp=10#NPC与HAC的酸性比例
k1=1.4
k4=0.4
k5=0
k6=0.00001



class Polymer():  #聚合物类(实体类)
    # 默认构造的链长度为0的聚合物（即是引发剂），支链为空，主链属性为True，如果主链属性为False则该链是某个链的支链
    def __init__(self,polymerLength=0,chainstatus=0,branch={},n_side=0,activate=True,n_branchnum=0):
        #取值为0 代表正常可以增长的链段，也就是活性链；
        #取值为1 代表是支链，
        #取值为2 代表是死链(死链死之前为主链) 其中死链和支链都不可以参与k4
        #取值为3 代表是回咬到了自己的任一深度的支链(死链死之前为主链) 其中死链和支链都不可以参与k4
        self.chainstatus = chainstatus  
        self.polymerLength = polymerLength  # 链长
        self.branch = {}  # 支链
        self.activate = activate # 是否是活性链，即是否有氨基
        self.n_side = 0 #初始的引发剂侧基数目为0
        self.n_branchnum = 0#初始支化点数目s为0

# 表征方法

#     def summary(): # 返回聚合物的基本信息
#         return mn,mw,pdi
#     def moleculearGPC(): # 将分子转化为GPC信号并绘制
#         transform
#          plot
    def structure(self):
        #根据分子数构造主链的结构
        self.polymerLength
#基元反应
'''
产生NCA单体，发生k1反应,
NPC和单体的数目减少，
随机选择一个聚合物，参加反应，并改变其结构

'''
import random

def acidification(): #酸化反应，没次反应循环都要发生一次
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    acidrate= (kp*n_NPC) / n_HA
    n_HA = kp*n_NPC + n_AcH; # 酸的总量
    k = 1.01
    a = k - 1
    b = k*n_HA - k*n_Polymer + 2*n_Polymer
    c = -n_Polymer**2
    n_activechain = (-b+(b**2-4*a*c)**0.5)/(2*a)
    n_HAremain = n_HA - n_Polymer + n_activechain
    activerate = n_activechain / n_Polymer
    n_NPCremain = n_HAremain * acidrate


def k1_reaction():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6

    # 修改数目
    n_NPC = n_NPC - 1
    n_NCA = n_NCA + 1
    n_phenol = n_phenol + 1  # 消耗的苯酚数目增加一个

'''
链增长
'''
def k4_reaction():  #
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    n_k4 = n_k4 + 1
    n_NCA = n_NCA - 1
    # 选择聚合物参加反应
    temp = choice(Polymer_list)
    # 判断聚合物的长度以及属性
    while temp.activate == False: #活性链才能够发生链增长反应
        # print('k4')
        temp = choice(Polymer_list)
    temp.polymerLength = temp.polymerLength + 1  #链长增加
    temp.n_side = temp.n_side + 1 # 参加反应的链的侧基数目+1
    n_side = n_side + 1 # 总的侧基数目+1

'''
链终止(变成死链或者支化)
'''

def k5_reaction():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    temp=choice(Polymer_list)
    while temp.polymerLength == 0:
        temp=choice(Polymer_list)#随机抽取一条链
    k = np.random.uniform(0,1)
    if k<=0.909:
        branch=get_all_branch(temp)
        if "temp.n_side" not in branch:
            temp.branch[temp.n_side]=temp#回咬最末端的侧基形成死链
    else:
        
        P = []
        branch = get_all_branch(temp)
        side=temp.n_side-1 #主链剩下侧基数
        for i in range(len(branch)):
            side=side+branch[i].n_side #加上侧链的侧基数
        print(side)
        if side > 0:
            for i in range(len(branch)):
                P.append(branch[i].n_side/side)
            P.append((temp.n_side-1)/side)
            
            temp3 = temp        
            accu_P = np.cumsum(P) #计算累积概率
            k = np.random.random()
            for i in range(len(accu_P)-1):
               if k <= accu_P[i]:
               #print(P[i])
                  temp3 = branch[i]
                  k=k+1
               #print(temp3)
            if temp3!=temp:      
               k = choice(sideloc(temp3))
            else:
               k=choice(sideloc2(temp3))
            temp3.branch[k] = temp
            temp3.n_side = temp3.n_side - 1
        else:
            temp.branch[temp.n_side]=temp
    temp.chainstatus = 2  #链变成死链
    temp.activate = False
    temp.n_side = temp.n_side - 1
    n_side = n_side - 1  
    n_Polymer = n_Polymer-1
    
    
def k6_reaction():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    Polymer_list2 = []
    Polymer_list2 = copy.deepcopy(Polymer_list)
    temp1 = choice(Polymer_list2)
    while temp1.activate == False: 
        temp1 = choice(Polymer_list2)
    # 然后根据各个链的侧基数目进行轮盘赌选择与它进行反应的链
    branch = get_all_branch(temp1)
    Polymer_list2.remove(temp1)
    for i in branch:
        Polymer_list2.remove(branch[i])
    P = []
    side=0
    for i in range(len(Polymer_list2)):
        side=side+Polymer_list2[i].n_side
    for i in range(len(Polymer_list2)):
            P.append((Polymer_list2[i].n_side)/side)
    accu_P = np.cumsum(P) #计算累积概率   
    k = np.random.random()
    for i in range(len(accu_P)):
        if k <= accu_P[i]:
            # print(P[i])
           temp2 = Polymer_list2[i]
           k = k + 1 # 跳出循环
    k = choice(sideloc(temp2)) 
    temp1.chainstatus = 1  # temp1变成支链不能增长
    temp2.branch[k] = temp1
    temp2.n_branchnum = temp2.n_branchnum + 1
    n_Polymer = n_Polymer-1
    temp2.n_side = temp2.n_side - 1
    temp1.activate=False
    n_side = n_side - 1
    
    
'''
获得一个链的所有支链对象（列表）
'''
def get_all_branch(Polymer):
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    branch = []
    for i in Polymer.branch:
        if Polymer.branch[i] == Polymer: #忽略占位的布尔变量
            pass
        else:
            branch.append(Polymer.branch[i])
            if len(Polymer.branch[i].branch) > 0:
                branch = branch + get_all_branch(Polymer.branch[i])
    return branch
'''
返回侧基的位置
'''


def sideloc(Polymer):
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    # 返回可以挂支链的位置
    loc = list(range(1,Polymer.polymerLength+1))
    for i in range(1,Polymer.polymerLength+1):
        if i in list(Polymer.branch.keys()):  # 挂有支链
            loc.remove(i)
    return loc


# 返回发生反应的函数名
def selectreaction():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    
    w_k1 = n_NPCanion * k1
    w_k4 = n_activechain * n_NCA * k4
    w_k5 = n_activechain * k5
    w_k6 = n_activechain * n_side * k6
    w = w_k1 + w_k4 + w_k5+w_k6
    p1 = w_k1 / w
    p4 = w_k4 / w
    p5 = w_k5 / w
    p6 = w_k6 / w
    k = np.random.random()#在0，1区间内产生随机数

    if k <= p1:
        return 'k1',w
    elif k <= p1 + p4:
        return 'k4',w
    elif k <= p1 + p4 + p5:
        return 'k5',w
    elif k <= p1 + p4 + p5 + p6:
        return 'k6',w
def alldead():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    for i in Polymer_list:
        if i.chainstatus == 0:
            return True
    return False 
# 计算聚合物的单体数和含有的支链数目
def calc_molecules(Polymer):
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6
    num = Polymer.polymerLength

    for i in Polymer.branch:
        if Polymer.branch[i] == Polymer or Polymer.branch[i] == False:
            pass
        else:
            num = num + calc_molecules(Polymer.branch[i])[0] #使用递归计算聚合物的单体数
    return num,len(Polymer.branch)+1
def calc_para():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6,n_NPC_init,n_Polymer_init,time
    w_single = 234
    w_initiator = 87
    M = 0
    M2 = 0
    Mw = 0
    Mn = 0
    pdi = 1
    n_chain = 0
    for i in Polymer_list:
        if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus == 3: #统计死链和活性链的个数
            if i.polymerLength > 1:
                n_chain = n_chain + 1
    # Mn = 每条链的质量/总链数 累加
    single_sum = 0
    for i in Polymer_list:

        if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus == 3:# 统计死链和活性链的个数
            if i.polymerLength >1:
                single_num,initiator_num = calc_molecules(i)
                M = M+ w_single*single_num
            # print('链上的单体数：%d'%single_num)
                single_sum = single_sum + single_num
                Mn = Mn + (w_single*single_num ) / n_chain # 数均分子量
                M2 = M2 + (w_single*single_num )**2# 重均分子量
    if M > 0:
        Mw = M2/M
        pdi = Mw/Mn
    

    # 
    return Mn,Mw,pdi,single_sum,n_chain


def calc_brach_num():#计算带有支链的分子数
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6,n_NPC_init,n_Polymer_init,time
    brach_num = 0
    n1=0
    m1=0
    for i in Polymer_list:
        if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus == 3: # 死链和活性链
            single_num,initiator_num = calc_molecules(i)
            if initiator_num>1:
                brach_num = brach_num + 1
                m1 = single_num*237
                n1= i.n_branchnum/single_num


    return brach_num



# 创建一个txt文件，文件名为当前时间,并将文件清空
filepath = r'/Users/viviantian/Downloads/'
now = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))
filename = filepath+now + 'result.txt'
f = open(filename, 'a')
#写入标题
title = 'Mn pdi n_branch_point n_brach_num_ratio n_NPC_ratio conversion time n_NCA_ratio n_NCA_ratio2 n_activechain_ratio n_ach_ratio n_brach_dead_ratio Mw n_Polymer_num_ratio,n_branchdegree'
f.write(title)


Mn_arr = []
Mw_arr = []
pdi_arr = []
NPC_arr = []
NCA_arr = []
n_Polymer_num_ratio_arr =[]
n_brach_num_ratio_arr = []
n_NPC_ratio_arr = []
n_NCA_ratio_arr = []
n_NCA_ratio2_arr = []
conversion_arr = []
n_activechain_ratio_arr = []
n_ach_ratio_arr = []
n_brach_dead_ratio_arr = []
n_branchdegree_arr = []
n_branchpoint_arr = []
w_molecular = [] # 分子量列表
Polymer_list = [Polymer() for i in range(n_Polymer)]  # 创建引发剂对象
conversion = 0
conversion_benchmark = 0
actualconnversion = [0,0.15854,0.27439,0.35976,0.51829,0.70122,0.82317,0.915]


while conversion<0.9 and alldead():
    acidification() # 酸化反应
    reaction,a0 = selectreaction()
    r = np.random.uniform(0, 1)# 在0，1区间内产生随机数
    # print(reaction)
    if reaction == 'k1':
        k1_reaction()
        T = T + 1/a0*np.log(1/r)
    elif reaction == 'k4':
        k4_reaction()
        T = T + 1/a0*np.log(1/r)
    elif reaction == 'k5':
        k5_reaction()
        T = T + 1/a0*np.log(1/r)
    elif reaction == 'k6':
        k6_reaction()
        T = T + 1/a0*np.log(1/r)
    #现有分子数/引发剂数目
    n_Polymer_num_ratio = n_Polymer/n_Polymer_init
    # 支化分子数目/现有分子数
    n_brach_num_ratio = calc_brach_num()/n_Polymer
    # 剩余NPC数目/初始NPC数目 –time
    n_NPC_ratio = n_NPC/n_NPC_init
    # 剩余NCA数目/初始NPC数目-time
    n_NCA_ratio = n_NCA/n_NPC_init
    # 反应掉的NCA数目/初始NPC数目—time====发生k4反应的次数/n_NPC_init
    n_NCA_ratio2 = n_k4/n_NPC_init
    # 活性链数目/现有分子数目—time or conversion
    n_activechain_ratio = n_activechain/n_Polymer
    # 酸化的链/现有分子数目—time or conversion
    n_ach_ratio = (n_Polymer-n_activechain)/n_Polymer
    # 支链+死链数目/现有分子数目-time or conversion
    n_brach_dead_ratio = (n_Polymer_init-n_Polymer)/ n_Polymer
    #支化度数目
    dpsum = 0
    for i in Polymer_list:

        if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus ==3:
            single_num,initiator_num = calc_molecules(i)
            if i.polymerLength != 0:
               dpsum = dpsum +i.n_branchnum/single_num
    n_branchdegree = dpsum/n_Polymer
    #支化点数目
    bpsum = 0
    for i in Polymer_list:
        if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus ==3:
              bpsum = bpsum +i.n_branchnum
    n_branchpoint = bpsum/n_Polymer


    n_brach_num_ratio_arr.append(n_brach_num_ratio)
    n_NCA_ratio_arr.append(n_NCA_ratio)
    n_NPC_ratio_arr.append(n_NPC_ratio)
    n_NCA_ratio2_arr.append(n_NCA_ratio2)
    n_activechain_ratio_arr.append(n_activechain_ratio)
    n_ach_ratio_arr.append(n_ach_ratio)
    n_brach_dead_ratio_arr.append(n_brach_dead_ratio)
    n_branchpoint_arr.append(n_branchpoint)
    n_branchdegree_arr.append(n_branchdegree)

    Mn,Mw,pdi,single_sum,n_chain = calc_para()
    Mn_arr.append(Mn)
    Mw_arr.append(Mw)
    pdi_arr.append(pdi)
    NPC_arr.append(n_NPC)
    NCA_arr.append(n_NCA)
    time_arr.append(T)
    w_molecular.append(single_sum)
    n_Polymer_num_ratio_arr.append(n_Polymer_num_ratio)
    
    conversion = n_phenol/n_NPC_init #反应速率
    # print(T)
    conversion_arr.append(conversion)
    if len(conversion_arr)== 1:
        #写入文件
        f.write('\n')
        f.writelines('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(Mn,pdi,n_branchpoint,n_brach_num_ratio,n_NPC_ratio,conversion,T,n_NCA_ratio,n_NCA_ratio2,n_activechain_ratio,n_ach_ratio,n_brach_dead_ratio,Mw,n_Polymer_num_ratio,n_branchdegree))
    else:
        if conversion - conversion_benchmark < 0.01:
            pass
        else:
            conversion_benchmark = conversion
            #写入文件
            f.write('\n')
            f.writelines('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(Mn,pdi,n_branchpoint,n_brach_num_ratio,n_NPC_ratio,conversion,T,n_NCA_ratio,n_NCA_ratio2,n_activechain_ratio,n_ach_ratio,n_brach_dead_ratio,Mw,n_Polymer_num_ratio,n_branchdegree))
    end = time.time()
    delta = end-start
#os.system('cls')
    # print('已经执行：%.2f%%\n大约剩余时间：%.2f分钟'%(((n_NPC_init-n_NPC)*100/n_NPC_init),((delta*n_NPC_init/(60*(n_NPC_init-n_NPC)))-delta/60)))
    print('已经执行：%.2f%%\n已经反应时间：%.2f分钟'%(((n_NPC_init-n_NPC)*100/n_NPC_init),(delta/60)))
f.close()

#gpc函数
filepath2 = r'/Users/viviantian/Downloads/'
now2 = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))
filename2 = filepath2 + now2 + 'gpc.txt'


def normfun(x,mu,sigma):
    pdf = np.exp(-((x-mu)**2)/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))
    return pdf
chainweight=[ ]

def plotgpc():
    global Polymer_list,n_Polymer,n_NPC,n_AcH,n_NCA,n_Ac,n_phenol,n_side,n_HA,n_activechain,n_NPCremain,k1,k4,k5,n_k4,n_k6,n_NPC_init,n_Polymer_init,time,Mn_arr,time_arr,s,inn,chainweight,w
#根据每个分子的Mn计算正态峰均值
    s=0
    inn=0
    n=0
    w = 0
    a = -2
    b = 22
    sigma = 0.5
    t = np.linspace(0,35,3000)
    gpc = [0 for i in t]
    for i in Polymer_list:
        if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus == 3:
            #统计死链和活性链的个数
            if i.polymerLength > 1:
                n = n + 1
    for i in Polymer_list:
        if i.chainstatus ==0 or i.chainstatus == 2 or i.chainstatus == 3:
            s, inn= calc_molecules(i)
            w = s*237/n
            chainweight.append(w)
    
    for i in chainweight:
        mu = -1*np.log10(i) + b
        #print(mu)cl
        gpc = gpc+ normfun(t,mu,sigma)

    return t, gpc

plt.figure(figsize=(10,10))
t,gpc = plotgpc()
    # update
f2 = open(filename2, 'a')
for j in gpc:
        #print(j)
    f2.write(str(j))
    f2.write('\n')
f2.close()
plt.plot(t,gpc)




print('n_NPC: %d\nn_NCA: %d'%(n_NPC,n_NCA))
print(alldead())
print('All have done!')

Mn,Mw,pdi,single_sum,n_chain = calc_para()
print('Mn: %.3f\nMw: %.3f\npdi: %.6f'%(Mn,Mw,pdi))
print('投入单体数: %d\n参加反应单体数: %d\n引发剂数目: %d'%(n_NPC_init,single_sum,n_Polymer_init))
print('链数：%d'%n_Polymer)

et=T/9
actualtime = [0,et,2*et,3*et,4.5*et,6*et,7.5*et,9*et]

filepath = r'/Users/viviantian/Downloads/'

filename = filepath+ 'branch.txt'
f = open(filename, 'a')
#写入标题
title = 'branchweight,branchdegree'
f.write(title)

m=0
n=0
single_num = 0
for i in Polymer_list:
    if i.chainstatus == 0 or i.chainstatus == 2 or i.chainstatus == 3: # 死链和活性链
        single_num,initiator_num = calc_molecules(i)
        if initiator_num>1:
            m=single_num*237
            n = i.n_branchnum
            f.write('\n')
            f.writelines('%s %s '%(m,n))

delta = end-start

# print(time_arr)
# 计算数目均分子量和重均分子量
# 主链的数目即为最后聚合物的数目
# plt.plot(time_arr,conversion_arr)
# plt.ylim(0,1)
# plt.xlabel('time')
# plt.ylabel('inversion rate')

plt.figure(figsize=(15,15))
plt.subplots_adjust(wspace =0.2, hspace =0.6)#调整子图间距
plt.subplot(331)
plt.plot(time_arr,Mn_arr)
plt.title('Mn-time')
plt.xlabel('time')
plt.ylabel('Mn')
plt.subplot(332)
plt.plot(conversion_arr,Mn_arr)
plt.title('Mn-conversion')
plt.xlabel('conversion')
plt.ylabel('Mn')
plt.subplot(333)
plt.plot(time_arr,pdi_arr)
plt.title('pdi-time')
plt.xlabel('time')
plt.ylabel('pdi')
plt.subplot(334)
plt.plot(conversion_arr,pdi_arr)
plt.title('pdi-conversion')
plt.xlabel('conversion')
plt.ylabel('pdi')
plt.subplot(335)
plt.plot(conversion_arr,n_branchpoint_arr)
plt.title('n_branchpoint-conversion')
plt.xlabel('conversion')
plt.ylabel('n_branchpoint')
plt.subplot(336)
plt.plot(time_arr,n_branchpoint_arr)
plt.title('n_branchpoint-time')
plt.xlabel('time')
plt.ylabel('n_branchpoint')

plt.subplot(337)
plt.plot(conversion_arr,n_brach_num_ratio_arr)
plt.title('n_brach_num_ratio-conversion')
plt.xlabel('conversion')
plt.ylabel('n_brach_num_ratio')

plt.subplot(338)
plt.plot(time_arr,n_NPC_ratio_arr)
plt.title('n_NPC_ratio-time')
plt.xlabel('time')
plt.ylabel('n_NPC_ratio')

plt.subplot(339)
plt.plot(time_arr,conversion_arr)
plt.ylim(0,1)
plt.xlabel('time')
plt.ylabel('conversion')
now = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))
plt.savefig('%s图1.png'%now,dpi=300)

plt.figure(figsize=(15,15))
plt.subplots_adjust(wspace =0.2, hspace =0.6)#调整子图间距
plt.subplot(331)
plt.plot(time_arr,n_NCA_ratio_arr)
plt.title('n_NCA_ratio-time')
plt.xlabel('time')
plt.ylabel('n_NCA_ratio')
plt.subplot(332)
plt.plot(time_arr,n_NCA_ratio2_arr)
plt.title('n_NCA_ratio2-time')
plt.xlabel('time')
plt.ylabel('n_NCA_ratio2')

plt.subplot(333)
plt.plot(time_arr,n_activechain_ratio_arr)
plt.title('n_activechain_ratio-time')
plt.xlabel('time')
plt.ylabel('n_activechain_ratio')
plt.subplot(334)
plt.plot(conversion_arr,n_activechain_ratio_arr)
plt.title('n_activechain_ratio-conversion')
plt.xlabel('conversion')
plt.ylabel('n_activechain_ratio')
plt.subplot(335)
plt.plot(time_arr,n_ach_ratio_arr)
plt.title('n_ach_ratio-time')
plt.xlabel('time')
plt.ylabel('n_ach_ratio')
plt.subplot(336)
plt.plot(conversion_arr,n_ach_ratio_arr)
plt.title('n_ach_ratio-conversion')
plt.xlabel('conversion')
plt.ylabel('n_ach_ratio')

plt.subplot(337)
plt.plot(time_arr,n_brach_dead_ratio_arr)
plt.title('n_brach_dead_ratio-time')
plt.xlabel('time')
plt.ylabel('n_brach_dead_ratio')
plt.subplot(338)
plt.plot(conversion_arr,n_brach_dead_ratio_arr)
plt.title('n_brach_dead_ratio-conversion')
plt.xlabel('conversion')
plt.ylabel('n_brach_dead_ratio')
plt.subplot(339)
plt.plot(time_arr,n_Polymer_num_ratio_arr)
plt.title('n_Polymer_num_ratio-time')
plt.xlabel('time')
plt.ylabel('n_Polymer_num_ratio')
now = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))
plt.savefig('%s图2.png'%now,dpi=300)

plt.figure(figsize=(15,15))
plt.subplots_adjust(wspace =0.2, hspace =0.6)
plt.subplot(331)
plt.plot(time_arr,n_branchdegree_arr)
plt.title('branchdegree-time')
plt.xlabel('time')
plt.ylabel('n_branchdegree')
plt.subplot(332)
plt.plot(conversion_arr,n_branchdegree_arr)
plt.title('branchdegree-conversion')
plt.xlabel('conversion')
plt.ylabel('n_branchdegree')
now = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time()))
plt.savefig('%s图3.png'%now,dpi=300)
plt.figure(figsize=(10,5))
y = actualconnversion
x = actualtime
plt.scatter(x,y)
plt.plot(time_arr,conversion_arr)
plt.ylim(0,1)
plt.xlabel('time')
plt.ylabel('conversion')
plt.title(conversion)
plt.show()
