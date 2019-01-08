# -*- coding: utf-8 -*-
import xlrd
import numpy
import numpy as np
import math
import numpy.linalg as LA
from sklearn.cluster import KMeans
import random
from sklearn import tree
a=open(r'data\Known disease-microbe association.xlsx') 
b=open(r'data\Symptom-based disease similarity.xlsx')
A=numpy.zeros((39,292))
SDM=numpy.zeros((39,39)) 
D1=numpy.zeros(910)   
D2=numpy.zeros(910)  
D3=numpy.zeros(910)  
D4=numpy.zeros(910)    
D5=numpy.zeros(910)  
D6=numpy.zeros(910)  
D7=numpy.zeros(910)  
D8=numpy.zeros(910)  
D9=numpy.zeros(910)  
D10=numpy.zeros(910)  
D11=numpy.zeros(910)  
D12=numpy.zeros(910)  
D13=numpy.zeros(910)  
D14=numpy.zeros(910)
D15=numpy.zeros(910)  
D16=numpy.zeros(910)  
D17=numpy.zeros(910)  
D18=numpy.zeros(910)   
D19=numpy.zeros(910)
D20=numpy.zeros(910)
D21=numpy.zeros(910)   
D22=numpy.zeros(910)  
D23=numpy.zeros(910)  
D24=numpy.zeros(910)    
D25=numpy.zeros(910)  
D26=numpy.zeros(910)  
D27=numpy.zeros(910)  
D28=numpy.zeros(910)  
D29=numpy.zeros(910)  
xlsx1=xlrd.open_workbook(r'data\Known disease-microbe association.xlsx')
sheet1=xlsx1.sheets()[0]
for i in range(450):
    s1=sheet1.row_values(i)                  
    m=int(s1[0])
    n=int(s1[1])
    A[m-1,n-1]=1                       #Adjacency matrix A
xlsx2=xlrd.open_workbook(r'data\symptom-based disease similarity.xlsx')
sheet2=xlsx2.sheets()[0]  
for i in range(61):
    s2=sheet2.row_values(i)                  
    z1=int(s2[0])
    z2=int(s2[1])
    z3=float(s2[2])
    SDM[z1-1,z2-1]=SDM[z2-1,z1-1]=z3
for i in range(39):
    SDM[i,i]=1                         #Symptom-based disease similarity SDM
#**************************************************************
C=np.asmatrix(A)
gamd=39/(LA.norm(C,'fro')**2);
kd=np.mat(np.zeros((39,39)))
km=np.mat(np.zeros((292,292)))
D=C*C.T;
for i in range(39):
        for j in range(i,39):
            kd[j,i]=np.exp(-gamd*(D[i,i]+D[j,j]-2*D[i,j]))
kd=kd+kd.T-np.diag(np.diag(kd))
KD=np.asarray(kd)                        #Gaussian interaction profile kernel similarity for disease
kd=[]                            
SD = (KD+SDM)/2
SD=np.asarray(SD)                        #Integrating symptom-based disease similarity SD
#***************************************************************
gamam = 292/(LA.norm(C,'fro')**2);
E=C.T*C;
for i in range(292):
    for j in range(i,292):
        km[i,j]=np.exp(-gamam*(E[i,i]+E[j,j]-2*E[i,j]))
km=km+km.T-np.diag(np.diag(km))
KM=np.asarray(km)                          #Gaussian interaction profile kernel similarity for microbe
km=[]                                    
#K_mean Clustering
unknown=[]                                 
known=[]                                  
for x in range(39):                        
    for y in range(292):
        if A[x,y]==0:                      
            unknown.append((x,y))
        else:
            known.append((x,y))             #Divide disease-microbe pairs into two parts based on whether they have known associations
major=[]
for z in range(10938):
    q=SD[unknown[z][0],:].tolist()+KM[unknown[z][1],:].tolist()
    major.append(q)
kmeans=KMeans(n_clusters=23, random_state=0).fit(major)
center=kmeans.cluster_centers_
center_x=[]
center_y=[]
for j in range(len(center)):
    center_x.append(center[j][0])
    center_y.append(center[j][1])
labels=kmeans.labels_
type1_x=[]
type1_y=[]
type2_x=[]
type2_y=[]
type3_x=[]
type3_y=[]
type4_x=[]
type4_y=[]
type5_x=[]
type5_y=[]
type6_x=[]
type6_y=[]
type7_x=[]
type7_y=[]
type8_x=[]
type8_y=[]
type9_x=[]
type9_y=[]
type10_x=[]
type10_y=[]
type11_x=[]
type11_y=[]
type12_x=[]
type12_y=[]
type13_x=[]
type13_y=[]
type14_x=[]
type14_y=[]
type15_x=[]
type15_y=[]
type16_x=[]
type16_y=[]
type17_x=[]
type17_y=[]
type18_x=[]
type18_y=[]
type19_x=[]
type19_y=[]
type20_x=[]
type20_y=[]
type21_x=[]
type21_y=[]
type22_x=[]
type22_y=[]
type23_x=[]
type23_y=[]
for i in range(len(labels)):
    if labels[i]==0:
        type1_x.append(unknown[i][0])
        type1_y.append(unknown[i][1])
    if labels[i]==1:
        type2_x.append(unknown[i][0])
        type2_y.append(unknown[i][1])
    if labels[i]==2:
        type3_x.append(unknown[i][0])
        type3_y.append(unknown[i][1])
    if labels[i]==3:
        type4_x.append(unknown[i][0])
        type4_y.append(unknown[i][1])
    if labels[i]==4:
        type5_x.append(unknown[i][0])
        type5_y.append(unknown[i][1])
    if labels[i]==5:
        type6_x.append(unknown[i][0])
        type6_y.append(unknown[i][1])
    if labels[i]==6:
        type7_x.append(unknown[i][0])
        type7_y.append(unknown[i][1])
    if labels[i]==7:
        type8_x.append(unknown[i][0])
        type8_y.append(unknown[i][1])
    if labels[i]==8:
        type9_x.append(unknown[i][0])
        type9_y.append(unknown[i][1])
    if labels[i]==9:
        type10_x.append(unknown[i][0])
        type10_y.append(unknown[i][1])
    if labels[i]==10:
        type11_x.append(unknown[i][0])
        type11_y.append(unknown[i][1])
    if labels[i]==11:
        type12_x.append(unknown[i][0])
        type12_y.append(unknown[i][1])
    if labels[i]==12:
        type13_x.append(unknown[i][0])
        type13_y.append(unknown[i][1])
    if labels[i]==13:
        type14_x.append(unknown[i][0])
        type14_y.append(unknown[i][1])
    if labels[i]==14:
        type15_x.append(unknown[i][0])
        type15_y.append(unknown[i][1])
    if labels[i]==15:
        type16_x.append(unknown[i][0])
        type16_y.append(unknown[i][1])
    if labels[i]==16:
        type17_x.append(unknown[i][0])
        type17_y.append(unknown[i][1])
    if labels[i]==17:
        type18_x.append(unknown[i][0])
        type18_y.append(unknown[i][1])
    if labels[i]==18:
        type19_x.append(unknown[i][0])
        type19_y.append(unknown[i][1])
    if labels[i]==19:
        type20_x.append(unknown[i][0])
        type20_y.append(unknown[i][1])
    if labels[i]==20:
        type21_x.append(unknown[i][0])
        type21_y.append(unknown[i][1])
    if labels[i]==21:
        type22_x.append(unknown[i][0])
        type22_y.append(unknown[i][1])
    if labels[i]==22:
        type23_x.append(unknown[i][0])
        type23_y.append(unknown[i][1])
type=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]                                       
mtype=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]                                     
dataSet=[]                                
for k1 in range(len(type1_x)):
    type[0].append((type1_x[k1],type1_y[k1]))      
for k2 in range(len(type2_x)):
    type[1].append((type2_x[k2],type2_y[k2]))
for k3 in range(len(type3_x)):
    type[2].append((type3_x[k3],type3_y[k3]))
for k4 in range(len(type4_x)):
    type[3].append((type4_x[k4],type4_y[k4]))
for k5 in range(len(type5_x)):
    type[4].append((type5_x[k5],type5_y[k5]))
for k6 in range(len(type6_x)):
    type[5].append((type6_x[k6],type6_y[k6]))
for k7 in range(len(type7_x)):
    type[6].append((type7_x[k7],type7_y[k7]))
for k8 in range(len(type8_x)):
    type[7].append((type8_x[k8],type8_y[k8]))     
for k9 in range(len(type9_x)):
    type[8].append((type9_x[k9],type9_y[k9]))
for k10 in range(len(type10_x)):
    type[9].append((type10_x[k10],type10_y[k10]))
for k11 in range(len(type11_x)):
    type[10].append((type11_x[k11],type11_y[k11]))
for k12 in range(len(type12_x)):
    type[11].append((type12_x[k12],type12_y[k12]))
for k13 in range(len(type13_x)):
    type[12].append((type13_x[k13],type13_y[k13]))
for k14 in range(len(type14_x)):
    type[13].append((type14_x[k14],type14_y[k14]))
for k15 in range(len(type15_x)):
    type[14].append((type15_x[k15],type15_y[k15]))     
for k16 in range(len(type16_x)):
    type[15].append((type16_x[k16],type16_y[k16]))
for k17 in range(len(type17_x)):
    type[16].append((type17_x[k17],type17_y[k17]))
for k18 in range(len(type18_x)):
    type[17].append((type18_x[k18],type18_y[k18]))
for k19 in range(len(type19_x)):
    type[18].append((type19_x[k19],type19_y[k19]))
for k20 in range(len(type20_x)):
    type[19].append((type20_x[k20],type20_y[k20]))
for k21 in range(len(type21_x)):
    type[20].append((type21_x[k21],type21_y[k21]))
for k22 in range(len(type22_x)):
    type[21].append((type22_x[k22],type22_y[k22]))
for k23 in range(len(type23_x)):
    type[22].append((type23_x[k23],type23_y[k23]))
for k in range(23):
    mtype[k]=random.sample(type[k],20)                 
for m2 in range(39):
    for n2 in range(292):
        for z2 in range(23):
            if (m2,n2) in mtype[z2]: 
                dataSet.append((m2,n2))                 
for m3 in range(39):
    for n3 in range(292):
        if A[m3,n3]==1:
            dataSet.append((m3,n3))                     #Building training samples                           
#Decision Tree
sumy1=numpy.zeros(910)
sumy2=numpy.zeros(910)
sumy3=numpy.zeros(910)
sumy4=numpy.zeros(910)
sumy5=numpy.zeros(910)
sumy6=numpy.zeros(910)
sumy7=numpy.zeros(910)
sumy8=numpy.zeros(910)
sumy9=numpy.zeros(910)
sumy10=numpy.zeros(910)
sumy11=numpy.zeros(910)
sumy12=numpy.zeros(910)
sumy13=numpy.zeros(910)
sumy14=numpy.zeros(910)
sumy15=numpy.zeros(910)
sumy16=numpy.zeros(910) 
sumy17=numpy.zeros(910)
sumy18=numpy.zeros(910)
sumy19=numpy.zeros(910)
sumy20=numpy.zeros(910)
sumy21=numpy.zeros(910)
sumy22=numpy.zeros(910)
sumy23=numpy.zeros(910)
sumy24=numpy.zeros(910)
sumy25=numpy.zeros(910)
sumy26=numpy.zeros(910)
sumy27=numpy.zeros(910)
sumy28=numpy.zeros(910)
sumy29=numpy.zeros(910)
sumy30=numpy.zeros(910)
x=[]                                                  
x1=[]                                                 
x2=[]
y=[]                                                                            
D=numpy.ones(910)*1.0/910.0                             #Initialization of training samples weights
for xx in dataSet:
    q=SD[xx[0],:].tolist()+KM[xx[1],:].tolist()
    x.append(q)                                         #Splicing of similarities of training samples 
    if (xx[0],xx[1]) in known:
        y.append(1)
    else:
        y.append(0)                                     #The labels of training samples
ys=numpy.array(y)
q=7
clf1=tree.DecisionTreeClassifier(max_depth=q)
clf1=clf1.fit(x,y,sample_weight=D)
v1=clf1.predict(x)                                      #Weak classifier predicts the labels of training samples
vs1=numpy.array(v1)
sumy1[vs1!=ys]=1                                  
sumd1=(sumy1*D).sum()
at1=math.log((1-sumd1)/sumd1)*0.5                       #The calculation of weights of weak classifiers in strong classifiers
z1=2*((sumd1*(1-sumd1))**0.5)       
for n3 in range(910):
    D1[n3]=D[n3]*numpy.exp(-at1*y[n3]*v1[n3])/z1        #The update of the weights of training samples 
clf2=tree.DecisionTreeClassifier(max_depth=q)
clf2=clf2.fit(x,y,sample_weight=D1)
v2=clf2.predict(x)                                       
vs2=numpy.array(v2)
sumy2[vs2!=ys]=1
sumd2=(sumy2*D1).sum()
at2=math.log((1-sumd2)/sumd2)*0.5         
z2=2*((sumd2*(1-sumd2))**0.5)       
for n3 in range(910):
    D2[n3]=D1[n3]*math.exp(-at2*y[n3]*v2[n3])/z2 
clf3=tree.DecisionTreeClassifier(max_depth=q)    
clf3=clf3.fit(x,y,sample_weight=D2)
v3=clf3.predict(x)              
vs3=numpy.array(v3)
sumy3[vs3!=ys]=1
sumd3=(sumy3*D2).sum()
at3=math.log((1-sumd3)/sumd3)*0.5         
z3=2*((sumd3*(1-sumd3))**0.5)       
for n3 in range(910):
    D3[n3]=D2[n3]*math.exp(-at3*y[n3]*v3[n3])/z3
clf4=tree.DecisionTreeClassifier(max_depth=q)
clf4=clf4.fit(x,y,sample_weight=D3)
v4=clf4.predict(x)                
vs4=numpy.array(v4)
sumy4[vs4!=ys]=1
sumd4=(sumy4*D3).sum()
at4=math.log((1-sumd4)/sumd4)*0.5        
z4=2*((sumd4*(1-sumd4))**0.5)       
for n3 in range(910):
    D4[n3]=D3[n3]*math.exp(-at4*y[n3]*v4[n3])/z4 
clf5=tree.DecisionTreeClassifier(max_depth=q)    
clf5=clf5.fit(x,y,sample_weight=D4)
v5=clf5.predict(x)               
vs5=numpy.array(v5)
sumy5[vs5!=ys]=1
sumd5=(sumy5*D4).sum()
at5=math.log((1-sumd5)/sumd5)*0.5        
z5=2*((sumd5*(1-sumd5))**0.5)       
for n3 in range(910):
    D5[n3]=D4[n3]*math.exp(-at5*y[n3]*v5[n3])/z5 
clf6=tree.DecisionTreeClassifier(max_depth=q)
clf6=clf6.fit(x,y,sample_weight=D5)
v6=clf6.predict(x)                
vs6=numpy.array(v6)
sumy6[vs6!=ys]=1
sumd6=(sumy6*D5).sum()
at6=math.log((1-sumd6)/sumd6)*0.5 
z6=2*((sumd6*(1-sumd6))**0.5)       
for n3 in range(910):
    D6[n3]=D5[n3]*math.exp(-at6*y[n3]*v6[n3])/z6 
clf7=tree.DecisionTreeClassifier(max_depth=q)
clf7=clf7.fit(x,y,sample_weight=D6)
v7=clf7.predict(x)               
vs7=numpy.array(v7)
sumy7[vs7!=ys]=1
sumd7=(sumy7*D6).sum()
at7=math.log((1-sumd7)/sumd7)*0.5        
z7=2*((sumd7*(1-sumd7))**0.5)       
for n3 in range(910):
    D7[n3]=D6[n3]*math.exp(-at7*y[n3]*v7[n3])/z7 
clf8=tree.DecisionTreeClassifier(max_depth=q)
clf8=clf8.fit(x,y,sample_weight=D7)
v8=clf8.predict(x)                
vs8=numpy.array(v8)
sumy8[vs8!=ys]=1
sumd8=(sumy8*D7).sum()
at8=math.log((1-sumd8)/sumd8)*0.5
z8=2*((sumd8*(1-sumd8))**0.5)
for n3 in range(910):
    D8[n3]=D7[n3]*math.exp(-at8*y[n3]*v8[n3])/z8
clf9=tree.DecisionTreeClassifier(max_depth=q)
clf9=clf9.fit(x,y,sample_weight=D8)
v9=clf9.predict(x)                                     
vs9=numpy.array(v9)
sumy9[vs9!=ys]=1
sumd9=(sumy9*D8).sum()
at9=math.log((1-sumd9)/sumd9)*0.5            
z9=2*((sumd9*(1-sumd9))**0.5)       
for n3 in range(910):
    D9[n3]=D8[n3]*math.exp(-at9*y[n3]*v9[n3])/z9           
clf10=tree.DecisionTreeClassifier(max_depth=q)
clf10=clf10.fit(x,y,sample_weight=D9)
v10=clf10.predict(x)                                       
vs10=numpy.array(v10)
sumy10[vs10!=ys]=1
sumd10=(sumy10*D9).sum()
at10=math.log((1-sumd10)/sumd10)*0.5         
z10=2*((sumd10*(1-sumd10))**0.5)       
for n3 in range(910):
    D10[n3]=D9[n3]*math.exp(-at10*y[n3]*v10[n3])/z10     
clf11=tree.DecisionTreeClassifier(max_depth=q)
clf11=clf11.fit(x,y,sample_weight=D10)
v11=clf11.predict(x)              
vs11=numpy.array(v11)
sumy11[vs11!=ys]=1
sumd11=(sumy11*D10).sum()
at11=math.log((1-sumd11)/sumd11)*0.5         
z11=2*((sumd11*(1-sumd11))**0.5)       
for n3 in range(910):
    D11[n3]=D10[n3]*math.exp(-at11*y[n3]*v11[n3])/z11 
clf12=tree.DecisionTreeClassifier(max_depth=q)
clf12=clf12.fit(x,y,sample_weight=D11)
v12=clf12.predict(x)                
vs12=numpy.array(v12)
sumy12[vs12!=ys]=1
sumd12=(sumy12*D11).sum()
at12=math.log((1-sumd12)/sumd12)*0.5        
z12=2*((sumd12*(1-sumd12))**0.5)       
for n3 in range(910):
    D12[n3]=D11[n3]*math.exp(-at12*y[n3]*v12[n3])/z12 
clf13=tree.DecisionTreeClassifier(max_depth=q)
clf13=clf13.fit(x,y,sample_weight=D12)
v13=clf13.predict(x)               
vs13=numpy.array(v13)
sumy13[vs13!=ys]=1
sumd13=(sumy13*D12).sum()
at13=math.log((1-sumd13)/sumd13)*0.5        
z13=2*((sumd13*(1-sumd13))**0.5)       
for n3 in range(910):
    D13[n3]=D12[n3]*math.exp(-at13*y[n3]*v13[n3])/z13    
clf14=tree.DecisionTreeClassifier(max_depth=q)
clf14=clf14.fit(x,y,sample_weight=D13)
v14=clf14.predict(x)                
vs14=numpy.array(v14)
sumy14[vs14!=ys]=1
sumd14=(sumy14*D13).sum()
at14=math.log((1-sumd14)/sumd14)*0.5 
z14=2*((sumd14*(1-sumd14))**0.5)       
for n3 in range(910):
    D14[n3]=D13[n3]*math.exp(-at14*y[n3]*v14[n3])/z14
clf15=tree.DecisionTreeClassifier(max_depth=q)
clf15=clf15.fit(x,y,sample_weight=D14)
v15=clf15.predict(x)               
vs15=numpy.array(v15)
sumy15[vs15!=ys]=1
sumd15=(sumy15*D14).sum()
at15=math.log((1-sumd15)/sumd15)*0.5        
z15=2*((sumd15*(1-sumd15))**0.5)       
for n3 in range(910):
    D15[n3]=D14[n3]*math.exp(-at15*y[n3]*v15[n3])/z15 
clf16=tree.DecisionTreeClassifier(max_depth=q)
clf16=clf16.fit(x,y,sample_weight=D15)
v16=clf16.predict(x)               
vs16=numpy.array(v16)
sumy16[vs16!=ys]=1
sumd16=(sumy16*D15).sum()
at16=math.log((1-sumd16)/sumd16)*0.5        
z16=2*((sumd16*(1-sumd16))**0.5)       
for n3 in range(910):
    D16[n3]=D15[n3]*math.exp(-at16*y[n3]*v16[n3])/z16
clf17=tree.DecisionTreeClassifier(max_depth=q)
clf17=clf17.fit(x,y,sample_weight=D16)
v17=clf17.predict(x)               
vs17=numpy.array(v17)
sumy17[vs17!=ys]=1
sumd17=(sumy17*D16).sum()
at17=math.log((1-sumd17)/sumd17)*0.5        
z17=2*((sumd17*(1-sumd17))**0.5)       
for n3 in range(910):
    D17[n3]=D16[n3]*math.exp(-at17*y[n3]*v17[n3])/z17
clf18=tree.DecisionTreeClassifier(max_depth=q)
clf18=clf18.fit(x,y,sample_weight=D17)
v18=clf18.predict(x)                
vs18=numpy.array(v18)
sumy18[vs18!=ys]=1
sumd18=(sumy18*D17).sum()
at18=math.log((1-sumd18)/sumd18)*0.5  
z18=2*((sumd18*(1-sumd18))**0.5)       
for n3 in range(910):
    D18[n3]=D17[n3]*math.exp(-at18*y[n3]*v18[n3])/z18 
clf19=tree.DecisionTreeClassifier(max_depth=q)
clf19=clf19.fit(x,y,sample_weight=D18)
v19=clf19.predict(x)                
vs19=numpy.array(v19)
sumy19[vs19!=ys]=1
sumd19=(sumy19*D18).sum()
at19=math.log((1-sumd19)/sumd19)*0.5  
z19=2*((sumd19*(1-sumd19))**0.5)       
for n3 in range(910):
    D19[n3]=D18[n3]*math.exp(-at19*y[n3]*v19[n3])/z19
clf20=tree.DecisionTreeClassifier(max_depth=q)
clf20=clf20.fit(x,y,sample_weight=D19)
v20=clf20.predict(x)                
vs20=numpy.array(v20)
sumy20[vs20!=ys]=1
sumd20=(sumy20*D19).sum()
at20=math.log((1-sumd20)/sumd20)*0.5                   
z20=2*((sumd20*(1-sumd20))**0.5)       
for n3 in range(910):
    D20[n3]=D19[n3]*math.exp(-at20*y[n3]*v20[n3])/z20
clf21=tree.DecisionTreeClassifier(max_depth=q)
clf21=clf21.fit(x,y,sample_weight=D20)
v21=clf21.predict(x)                                      
vs21=numpy.array(v21)
sumy21[vs21!=ys]=1                                  
sumd21=(sumy21*D20).sum()
at21=math.log((1-sumd21)/sumd21)*0.5                       
z21=2*((sumd21*(1-sumd21))**0.5)       
for n3 in range(910):
    D21[n3]=D20[n3]*numpy.exp(-at21*y[n3]*v21[n3])/z21      
clf22=tree.DecisionTreeClassifier(max_depth=q)
clf22=clf22.fit(x,y,sample_weight=D21)
v22=clf22.predict(x)                                       
vs22=numpy.array(v22)
sumy22[vs22!=ys]=1
sumd22=(sumy22*D21).sum()
at22=math.log((1-sumd22)/sumd22)*0.5         
z22=2*((sumd22*(1-sumd22))**0.5)       
for n3 in range(910):
    D22[n3]=D21[n3]*math.exp(-at22*y[n3]*v22[n3])/z22 
clf23=tree.DecisionTreeClassifier(max_depth=q)    
clf23=clf23.fit(x,y,sample_weight=D22)
v23=clf23.predict(x)              
vs23=numpy.array(v23)
sumy23[vs23!=ys]=1
sumd23=(sumy23*D22).sum()
at23=math.log((1-sumd23)/sumd23)*0.5         
z23=2*((sumd23*(1-sumd23))**0.5)       
for n3 in range(910):
    D23[n3]=D22[n3]*math.exp(-at23*y[n3]*v23[n3])/z23
clf24=tree.DecisionTreeClassifier(max_depth=q)
clf24=clf24.fit(x,y,sample_weight=D23)
v24=clf24.predict(x)                
vs24=numpy.array(v24)
sumy24[vs24!=ys]=1
sumd24=(sumy24*D23).sum()
at24=math.log((1-sumd24)/sumd24)*0.5        
z24=2*((sumd24*(1-sumd24))**0.5)       
for n3 in range(910):
    D24[n3]=D23[n3]*math.exp(-at24*y[n3]*v24[n3])/z24 
clf25=tree.DecisionTreeClassifier(max_depth=q)    
clf25=clf25.fit(x,y,sample_weight=D24)
v25=clf25.predict(x)               
vs25=numpy.array(v25)
sumy25[vs25!=ys]=1
sumd25=(sumy25*D24).sum()
at25=math.log((1-sumd25)/sumd25)*0.5        
z25=2*((sumd25*(1-sumd25))**0.5)       
for n3 in range(910):
    D25[n3]=D24[n3]*math.exp(-at25*y[n3]*v25[n3])/z25 
clf26=tree.DecisionTreeClassifier(max_depth=q)
clf26=clf26.fit(x,y,sample_weight=D25)
v26=clf26.predict(x)                
vs26=numpy.array(v26)
sumy26[vs26!=ys]=1
sumd26=(sumy26*D25).sum()
at26=math.log((1-sumd26)/sumd26)*0.5 
z26=2*((sumd26*(1-sumd26))**0.5)       
for n3 in range(910):
    D26[n3]=D25[n3]*math.exp(-at26*y[n3]*v26[n3])/z26 
clf27=tree.DecisionTreeClassifier(max_depth=q)
clf27=clf27.fit(x,y,sample_weight=D26)
v27=clf27.predict(x)               
vs27=numpy.array(v27)
sumy27[vs27!=ys]=1
sumd27=(sumy27*D26).sum()
at27=math.log((1-sumd27)/sumd27)*0.5        
z27=2*((sumd27*(1-sumd27))**0.5)       
for n3 in range(910):
    D27[n3]=D26[n3]*math.exp(-at27*y[n3]*v27[n3])/z27 
clf28=tree.DecisionTreeClassifier(max_depth=q)
clf28=clf28.fit(x,y,sample_weight=D27)
v28=clf28.predict(x)                
vs28=numpy.array(v28)
sumy28[vs28!=ys]=1
sumd28=(sumy28*D27).sum()
at28=math.log((1-sumd28)/sumd28)*0.5
z28=2*((sumd28*(1-sumd28))**0.5)
for n3 in range(910):
    D28[n3]=D27[n3]*math.exp(-at28*y[n3]*v28[n3])/z28
clf29=tree.DecisionTreeClassifier(max_depth=q)
clf29=clf29.fit(x,y,sample_weight=D28)
v29=clf29.predict(x)                                     
vs29=numpy.array(v29)
sumy29[vs29!=ys]=1
sumd29=(sumy29*D28).sum()
at29=math.log((1-sumd29)/sumd29)*0.5            
z29=2*((sumd29*(1-sumd29))**0.5)       
for n3 in range(910):
    D29[n3]=D28[n3]*math.exp(-at29*y[n3]*v29[n3])/z29           
clf30=tree.DecisionTreeClassifier(max_depth=q)
clf30=clf30.fit(x,y,sample_weight=D29)
v30=clf30.predict(x)                                       
vs30=numpy.array(v30)
sumy30[vs30!=ys]=1
sumd30=(sumy30*D29).sum()
at30=math.log((1-sumd30)/sumd30)*0.5                         #Training process is completed
for yy in unknown:
    q1=SD[yy[0],:].tolist()+KM[yy[1],:].tolist()           
    x1.append(q1)  
fs=clf1.predict_proba(x1)*at1+clf2.predict_proba(x1)*at2+clf3.predict_proba(x1)*at3+clf4.predict_proba(x1)*at4+clf5.predict_proba(x1)*at5+clf6.predict_proba(x1)*at6+clf7.predict_proba(x1)*at7+clf8.predict_proba(x1)*at8+clf9.predict_proba(x1)*at9+clf10.predict_proba(x1)*at10+clf11.predict_proba(x1)*at11+clf12.predict_proba(x1)*at12+clf13.predict_proba(x1)*at13+clf14.predict_proba(x1)*at14+clf15.predict_proba(x1)*at15+clf16.predict_proba(x1)*at16+clf17.predict_proba(x1)*at17+clf18.predict_proba(x1)*at18+clf19.predict_proba(x1)*at19+clf20.predict_proba(x1)*at20+clf21.predict_proba(x1)*at21+clf22.predict_proba(x1)*at22+clf23.predict_proba(x1)*at23+clf24.predict_proba(x1)*at24+clf25.predict_proba(x1)*at25+clf26.predict_proba(x1)*at26+clf27.predict_proba(x1)*at27+clf28.predict_proba(x1)*at28+clf29.predict_proba(x1)*at29+clf30.predict_proba(x1)*at30
px1=fs[:,1].tolist()                                         #The scores of unknown samples                       
#print(px1)
sy=[0]*10938                                                       
xlsx7=xlrd.open_workbook(r'data\Disease number.xlsx')
xlsx8=xlrd.open_workbook(r'data\Microbe number.xlsx')
sheet7=xlsx7.sheets()[0]                                 
sheet8=xlsx8.sheets()[0]  
px1=numpy.matrix(px1)
Sampleranking=numpy.argsort(-px1).tolist()    
Sampleranking=Sampleranking[0]  
f=open('Prediction results for all unknown samples.txt','a+')
f.writelines(['disease','\t','microbe','\t','Score','\n'])
f.close()
for i in range(10938):                     
    a=fs[:,1][Sampleranking[i]]
    s7=sheet7.row_values(unknown[Sampleranking[i]][0])                   
    s8=sheet8.row_values(unknown[Sampleranking[i]][1])                    
    f=open('Prediction results for all unknown samples.txt','a+')
    f.writelines([s7[1],'\t',s8[1],'\t',str(a),'\n'])
    f.close()                                                #Getting the prediction results for all unknown samples  

   
