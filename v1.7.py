import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from obspy import Trace
from pykalman import KalmanFilter


'''
conda install -c obspy obspy 
conda install -c conda-forge matplotlib 
'''
#Read data from .csv given by
def read():
    #f=open('./imu_log/27cm/30sPlus27cm3times01.csv')
    #f=open('./imu_log/29cm/X_29cm_3times_02.csv')
    f=open('./imu_log/14.5cm/14.5X2X02.csv')
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
    print()

def correctOrientation(radians): #not finished
    x_R = np.asarray([[1, 0, 0],
                   [0, math.cos(radians[0]), -1*math.sin(radians[0])],
                   [0, math.sin(radians[0]), math.cos(radians[0])]], dtype=np.float64)
    y_R = np.asarray([[math.cos(radians[1]), 0, math.sin(radians[1])],
                   [0, 1, 0],
                   [-1*math.sin(radians[1]), 0, math.cos(radians[1])]], dtype=np.float64)
    z_R = np.asarray([[math.cos(radians[2]), -1 * math.sin(radians[2]), 0],
                   [math.sin(radians[2]), math.cos(radians[2]), 0],
                   [0, 0, 1]], dtype=np.float64)
    Rotation = np.dot( np.dot(x_R, y_R), z_R)
    return Rotation

def calVelocityAndDisplacement(anvNacc, sampleRate,vel_record, dsp_record):
    for i in range(1, len(anvNacc)+1):
        period = 1/sampleRate
        dsp_record[i][0] = vel_record[i-1][0]*period+(1/2)*anvNacc[i-1][3]*(period**2)
        dsp_record[i][1] = vel_record[i-1][1]*period+(1/2)*anvNacc[i-1][4]*(period**2)
        dsp_record[i][2] = vel_record[i-1][2]*period+(1/2)*anvNacc[i-1][5]*(period**2)
        vel_record[i][0] = vel_record[i-1][0] + anvNacc[i-1][3]*period
        vel_record[i][1] = vel_record[i-1][1] + anvNacc[i-1][4]*period
        vel_record[i][2] = vel_record[i-1][2] + anvNacc[i-1][5]*period
    for i in range(1, len(anvNacc)+1):
        dsp_record[i] += dsp_record[i-1]

def showPlot_1_Axes(array, name, unit, lowerbound=0, upperbound=0): #array with 1 parameters
    plt.figure()
    plt.title(name)
    plt.plot(array,color='red',label='rec')
    plt.legend()

    plt.xlabel("(ms)")
    plt.ylabel('('+unit+')')
    if upperbound != 0 and lowerbound != 0:
        plt.ylim(lowerbound, upperbound)
    plt.show()


def correctCoordinate(anvNacc, sample_num):
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


def calByMatrix(acc, sampleRate):
    dt = 1/sampleRate
    F = [[1, dt, 0.5*dt**2], 
     [0,  1,       dt],
     [0,  0,        1]]
    F = np.asarray(F, dtype=np.float64)
    X = np.zeros([3, 1],dtype=np.float64)
    rec = np.zeros([len(acc)+1, 3], dtype=np.float64)
    for i in range(1, len(acc)+1):
        X[2] = acc[i-1]
        X = np.dot(F, X)
        rec[i] = np.swapaxes(X,0,1)
    print(X[0])
    return rec


def main():
    '''
    anvNacc = read() #data input
    aver_sample_num = 30000 #1000 sample as 1 sec.
    #anvNacc = correctCoordinate(anvNacc, aver_sample_num)
    anvNacc = anvNacc[:30000]
    sns.set()
    #x = np.random.randn(100)
    #x = np.asarray([1,1,1,3,5,6,6])
    
    sns.distplot(anvNacc[:,3])
    plt.show()
    '''
    sampleNumber = 18000
    sampleRate = 1
    anvNacc = np.zeros([sampleNumber,6],dtype=np.float64)
    rad = np.linspace(0, 6*np.pi, sampleNumber)
    coswave = np.cos(rad)
    #app = np.asarray([0,1,2,3,4,5,4,3,2,1,0,-1,-2,-3,-4,-5,-4,-3,-2,-1,0])
    anvNacc[:,3] += coswave

    ########## Filter Field ##########
    '''
    tr=Trace()
    tr.stats['sampling_rate']=1000
    for i in range(6):
        tr.data=anvNacc[:,i]
        #tr.filter("bandpass", freqmin =0.027,freqmax 100)
        tr.filter("highpass", freq =0.1)
        #tr.filter("lowpass", freq=100)
        anvNacc[:,i] = tr.data
    '''
    ########## Block of Computing Velocity & Displacement ##########
    '''
    record = calByMatrix(anvNacc[:,3], 1000)
    showPlot_1_Axes(anvNacc[:,3], "Acceleration", "m/s2")
    showPlot_1_Axes(record[:,1], "Velocity", "m/s")
    showPlot_1_Axes(record[:,0], "Displacement", "m")
    print(record[999, 2])
    print(record[999, 1])
    print(record[999, 0])
    '''
    '''
    anvNacc = np.array(
            [[0,0,0,2,3,2],
            [0,0,0,4,4,1],
            [0,0,0,5,1,7],
            [0,0,0,3,0,5]],dtype=np.float64)
    '''
    vel_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the velocity
    dsp_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the displacement
    calVelocityAndDisplacement(anvNacc, sampleRate, vel_record, dsp_record)

    showPlot_1_Axes(anvNacc[:,3], "Acceleration", "m/s2")
    showPlot_1_Axes(vel_record[:,0], "Velocity", "m/s")
    showPlot_1_Axes(dsp_record[:,0], "Displacement", "m")
    print(anvNacc[len(anvNacc)-1, 3])
    print(vel_record[len(anvNacc)])
    print(dsp_record[len(anvNacc)])



    
    mu, variance = 0, 1.5 # mean and standard deviation
    s = np.random.normal(mu, variance, sampleNumber) #Create Gaussian Noise
    anvNacc[:,3] += s

    vel_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the velocity
    dsp_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the displacement
    calVelocityAndDisplacement(anvNacc, sampleRate, vel_record, dsp_record)

    showPlot_1_Axes(anvNacc[:,3], "Acceleration", "m/s2")
    showPlot_1_Axes(vel_record[:,0], "Velocity", "m/s")
    showPlot_1_Axes(dsp_record[:,0], "Displacement", "m")
    print(anvNacc[len(anvNacc)-1, 3])
    print(vel_record[len(anvNacc)])
    print(dsp_record[len(anvNacc)])

    tr=Trace()
    tr.stats['sampling_rate']=sampleRate
    tr.data=anvNacc[:,3]
    tr.filter("bandpass", freqmin =0.08,freqmax = 100)
    anvNacc[:,3] = tr.data

    vel_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the velocity
    dsp_record = np.zeros([len(anvNacc)+1, 3],dtype=np.float64) # To record the displacement
    calVelocityAndDisplacement(anvNacc, sampleRate, vel_record, dsp_record)

    showPlot_1_Axes(anvNacc[:,3], "Acceleration", "m/s2")
    showPlot_1_Axes(vel_record[:,0], "Velocity", "m/s")
    showPlot_1_Axes(dsp_record[:,0], "Displacement", "m")
    print(anvNacc[len(anvNacc)-1, 3])
    print(vel_record[len(anvNacc)])
    print(dsp_record[len(anvNacc)])
    #Plot Gaussian Noise

    #plt.xlim(0,1000)
    #plt.ylim(-20,20)
    plt.xlabel("(ms)")
    plt.ylabel("(m/s^2)")
    plt.plot(np.linspace(1, sampleNumber, num = sampleNumber),s )
    plt.show()
    
    #Plot Gaussian Distribution

    count, bins, ignored = plt.hist(s, 300, normed=True)
    plt.plot(bins, 1/(variance * np.sqrt(2 * np.pi)) * 
            np.exp( - (bins - mu)**2 / (2 * variance**2) ), 
            linewidth=2, color='r')
    plt.xlabel("(m/s^2)")
    plt.ylabel("(proportion)")    
    plt.show()
    

if __name__ == '__main__':
    main()