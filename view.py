import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from obspy import Trace


def read():

    f=open('./imu_log/29cm/X_29cm_3times_02.csv')
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


plt.figure()

plt.title('Raw Data of Acceleration')
plt.plot(anvNacc[:,3],color='green',label='x_acc')
plt.plot(anvNacc[:,4],color='red',label='y_acc')
plt.plot(anvNacc[:,5],color='blue',label='z_acc')

plt.legend()

plt.xlabel('ms')
plt.ylabel('m/s2')
#plt.ylim(.8, -.8)
#plt.ylim(-9, -10.6)
plt.show()



########## Filter Field ##########
tr=Trace()
tr.stats['sampling_rate']=1000
for i in range(6):
    tr.data=anvNacc[:,i]
    tr.filter("bandpass", freqmin =0.1,freqmax =100) 
    #tr.filter("highpass", freq =0.1) 
    #tr.filter("lowpass", freq=100)
    anvNacc[:,i] = tr.data
########## End of Filter Field ##########
	
plt.figure()

plt.title('Pruned Data of Acceleration')
plt.plot(anvNacc[:,3],color='green',label='x_acc')
plt.plot(anvNacc[:,4],color='red',label='y_acc')
plt.plot(anvNacc[:,5],color='blue',label='z_acc')

plt.legend()

plt.xlabel('ms')
plt.ylabel('m/s2')
#plt.ylim(-1, 1)
plt.show()

