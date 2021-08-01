import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from obspy import Trace

def read():

    f=open('noFilter15cm3times.csv')
    lines=f.readlines()[6:]

    data=list()
    for line in lines:
        temp=np.array(line.strip().split(','))
        if temp[-1]=='OK':
            temp2=temp[:-6].astype(np.float64)
            temp2=np.delete(temp2,1,0)
            temp2=np.delete(temp2,0,0)
            data.append(temp2)

    data = np.asarray(data)
    f.close()
    return data

anvNacc = read() #data input
anvNacc[:,5]+=9.8
#anvNacc[:,5]=0

tr=Trace()
tr.stats['sampling_rate']=1000
for i in range(6):

    
    tr.data=anvNacc[:,i]
    tr.filter("bandpass", freqmin =0.1,freqmax =100) 
    #tr.filter("lowpass", freq=0.1) 
    anvNacc[:,i] = tr.data


cur_an = np.zeros(3, dtype=np.float64)
X = np.zeros(6, dtype=np.float64).reshape(6,1)
A = np.asarray([[1, 0, 0, 0.001, 0, 0],
               [0, 1, 0, 0, 0.001, 0],
               [0, 0, 1, 0, 0, 0.001],
               [0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 1]], dtype=np.float64)
Orientation = np.eye(3, dtype=np.float64)
U = np.zeros(3, dtype=np.float64).reshape(3, 1)
B= np.vstack((np.dot(np.eye(3, dtype=np.float64), 0.001*0.001/2), np.dot(np.eye(3, dtype=np.float64), 0.001)))
X_record = np.array([0,0,0], dtype=np.float64)

cur_ans = list()
for i in range(len(anvNacc)):

    if i == 0:
        cur_an = (np.zeros(3, dtype=np.float64) + anvNacc[i][0:3])/2*0.001
    else:
        cur_an = (anvNacc[i-1][0:3] + anvNacc[i][0:3])/2*0.001

    cur_ans.append(cur_an)
cur_ans=np.array(cur_ans)
'''for i in range(3):

    
    tr.data=cur_ans[:,i]
    tr.filter("bandpass", freqmin =0.1,freqmax =100) 
    #tr.filter("lowpass", freq=0.1) 
    cur_ans[:,i] = tr.data'''

xs=list()
Rotations=list()
for i in range(len(anvNacc)):
    if i == 0:
        cur_an = (np.zeros(3, dtype=np.float64) + anvNacc[i][0:3])/2*0.001
    else:
        cur_an = (anvNacc[i-1][0:3] + anvNacc[i][0:3])/2*0.001
    #cur_an=cur_ans[i]
    u_R = np.asarray([[1, 0, 0],
                   [0, math.cos(cur_an[0]), -1*math.sin(cur_an[0])],
                   [0, math.sin(cur_an[0]), math.cos(cur_an[0])]], dtype=np.float64)
    v_R = np.asarray([[math.cos(cur_an[1]), 0, math.sin(cur_an[1])],
                   [0, 1, 0],
                   [-1*math.sin(cur_an[1]), 0, math.cos(cur_an[1])]], dtype=np.float64)
    w_R = np.asarray([[math.cos(cur_an[2]), -1 * math.sin(cur_an[2]), 0],
                   [math.sin(cur_an[2]), math.cos(cur_an[2]), 0],
                   [0, 0, 1]], dtype=np.float64)
    
    
    Rotation = np.dot(w_R, np.dot(v_R,u_R))
    Rotations.append(Rotation)
    '''
    try:
        Rotation = np.linalg.inv(Rotation)
    except np.linalg.LinAlgError:
        # Not invertible. Skip this one.
        pass
    '''
    #X = np.dot(A, X) + np.dot( B, 0.5*( ( np.dot( np.dot(Orientation, Rotation), anvNacc[i][3:6].reshape(3, 1))) +  np.dot(Orientation, U)) )
    #Rotation = np.eye(3)
    #Orientation = np.eye(3)
    #X = np.dot(A, X) + np.dot( B, 0.5*( ( np.dot( np.dot(Orientation, Rotation), anvNacc[i][3:6].reshape(3, 1))) +  np.dot(Orientation, U)) )
    #x=0.5*( ( np.dot( np.dot(Rotation,Orientation), anvNacc[i][3:6].reshape(3, 1))) +  np.dot(Orientation, U))
    Orientation = np.dot( Rotation,Orientation)
    x=0.5*( np.dot(Orientation,(anvNacc[i][3:6].reshape(3, 1)+U)))
    #print (Orientation.shape)
    xs.append(np.squeeze(x))
    #print (np.squeeze(x).shape)
    #exit()
    #X_record = np.vstack((X_record,X[0:3].reshape(1,3)))
    
    U = anvNacc[i][3:6].reshape(3, 1)

xs=np.array(xs)
#print (xs.shape)
'''for i in range(3):

    tr.data=xs[:,i]
    tr.filter("highpass", freq =0.1) 
    #tr.filter("lowpass", freq=0.1) 
    xs[:,i] = tr.data
'''
#print (xs.shape)
#exit()
X_record = np.array([0,0,0], dtype=np.float64)
X1=np.zeros(3)
for i in range(len(xs)):

    Rot=np.block([
        [np.eye(3),np.zeros((3,3))],
        [np.zeros((3,3)),Rotations[i]]
    ])
    X = np.dot(A, np.dot(Rot,X)) + np.dot( B, xs[i].reshape(3,1) )
    '''if i==0:
        X1+=(xs[i]/2)*0.001
    else:
        X1+=((xs[i]+xs[i-1])/2)*0.001'''
    
    #print (X.shape)
    #exit()
    X_record = np.vstack((X_record,X[0:3].reshape(1,3)))

'''
ax = plt.subplot(111, projection='3d')
ax.scatter(X_record[:,0], X_record[:,1], X_record[:,2], c='y')

print(X)
ax.set_zlabel('Z')  
ax.set_ylabel('Y')
ax.set_xlabel('X')
plt.show()


plt.figure()
plt.plot(X_record[:,2])
plt.show()
'''
plt.figure()

plt.plot(X_record[:,0],X_record[:,1])

plt.show()
