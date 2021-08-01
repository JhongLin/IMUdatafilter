import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from obspy import Trace

'''
conda install -c obspy obspy 
conda install -c conda-forge matplotlib 


'''
#Read data from .csv given by
def read():
    f=open('./imu_log/27cm/30sPlus27cm3times01.csv')
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
#Utilize the gravity vector to correct 3-axes component
def initialAngleCorrection(vx, vy, vz):
    Gravity = -9.80665
    radians = np.zeros(3, dtype=np.float64)
    radians[0] = math.acos(vx/Gravity)
    radians[1] = math.acos(vy/Gravity)
    radians[2] = math.acos(vz/Gravity)
    return radians

def gravityOrientation(radians):
    u_R = np.asarray([[1, 0, 0],
                   [0, math.cos(radians[0]), -1*math.sin(radians[0])],
                   [0, math.sin(radians[0]), math.cos(radians[0])]], dtype=np.float64)
    v_R = np.asarray([[math.cos(radians[1]), 0, math.sin(radians[1])],
                   [0, 1, 0],
                   [-1*math.sin(radians[1]), 0, math.cos(radians[1])]], dtype=np.float64)
    w_R = np.asarray([[math.cos(radians[2]), -1 * math.sin(radians[2]), 0],
                   [math.sin(radians[2]), math.cos(radians[2]), 0],
                   [0, 0, 1]], dtype=np.float64)
    Rotation = np.dot( np.dot(u_R, v_R), w_R)
    return Rotation

def velocity(anvNacc, sample_num):
    v = np.zeros([len(anvNacc)+1, 3], dtype=np.float64)
    for i in range(1, len(anvNacc)+1):
        v[i][0] = v[i-1][0]+anvNacc[i-1][3]*(0.001)
        v[i][1] = v[i-1][1]+anvNacc[i-1][4]*(0.001)
        v[i][2] = v[i-1][2]+anvNacc[i-1][5]*(0.001)
    return v

def displacment():
    print()

def correctCoordinate(anvNacc, sample_num):
    Gravity = -9.80665

    if sample_num == 0:
        return anvNacc
    
    sum_xyz = np.zeros(3,dtype=np.float64)
    aver_xyz = np.zeros(3,dtype=np.float64)
    for i in range(sample_num):
        sum_xyz[0] += anvNacc[i][3]
        sum_xyz[1] += anvNacc[i][4]
        sum_xyz[2] += anvNacc[i][5]

    aver_xyz[0] = sum_xyz[0] / sample_num
    aver_xyz[1] = sum_xyz[1] / sample_num
    aver_xyz[2] = sum_xyz[2] / sample_num
    #rads = initialAngleCorrection(x_aver, y_aver, z_aver)
    anvNacc[:,3] -= aver_xyz[0]
    anvNacc[:,4] -= aver_xyz[1]
    anvNacc[:,5] -= aver_xyz[2]
    return anvNacc[sample_num:]
def main():
    '''
    anvNacc = read() #data input
    aver_sample_num = 30000 #1000 sample as 1 sec.
    anvNacc = correctCoordinate(anvNacc, aver_sample_num)
    '''
    '''
    ########## Filter Field ##########
    tr=Trace()
    tr.stats['sampling_rate']=1000
    for i in range(6):
        tr.data=anvNacc[:,i]
        #tr.filter("bandpass", freqmin =0.1,freqmax =100) 
        tr.filter("highpass", freq =0.1) 
        #tr.filter("lowpass", freq=0.1)
        anvNacc[:,i] = tr.data
    ########## Filter Field ##########
    '''
    anvNacc = np.zeros([18000,6],dtype=np.float64)
    rad = np.linspace(0, 2*np.pi, 18000)
    coswave = np.cos(rad)
    anvNacc[:,3] = coswave

    vel_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the velocity
    dsp_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the displacement
    
    for i in range(1, len(anvNacc)+1):
        dsp_record[i][0] = vel_record[i-1][0]*(0.001)+(1/2)*anvNacc[i-1][3]*(0.000001)
        dsp_record[i][1] = vel_record[i-1][1]*(0.001)+(1/2)*anvNacc[i-1][4]*(0.000001)
        dsp_record[i][2] = vel_record[i-1][2]*(0.001)+(1/2)*anvNacc[i-1][5]*(0.000001)
        vel_record[i][0] = vel_record[i-1][0] + anvNacc[i-1][3]*(0.001)
        vel_record[i][1] = vel_record[i-1][1] + anvNacc[i-1][4]*(0.001)
        vel_record[i][2] = vel_record[i-1][2] + anvNacc[i-1][5]*(0.001)
        #print('%f %f %f' %(vel_record[i][0], vel_record[i][1], vel_record[i][2]))
    print('----------------------------')
    #print(dsp_record[0])
    for i in range(1, len(anvNacc)+1):
        dsp_record[i] += dsp_record[i-1]
        #print(dsp_record[i])
    
    print(dsp_record[len(anvNacc)])
    
    plt.figure()
    plt.title('Displacement')
    plt.plot(dsp_record[:,0],color='green',label='x_rec')
    plt.plot(dsp_record[:,1],color='red',label='y_rec')
    plt.plot(dsp_record[:,2],color='blue',label='z_rec')

    plt.legend()

    plt.xlabel('(ms)')
    plt.ylabel('(m)')
    #plt.ylim(0, 10)
    #plt.ylim(-9, -10.6)
    plt.show()
    
if __name__ == '__main__':
    main()