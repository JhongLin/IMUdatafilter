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

def calVelocityAndDisplacement(anvNacc, vel_record, dsp_record):
    for i in range(1, len(anvNacc)+1):
        dsp_record[i][0] = vel_record[i-1][0]*(0.001)+(1/2)*anvNacc[i-1][3]*(0.000001)
        dsp_record[i][1] = vel_record[i-1][1]*(0.001)+(1/2)*anvNacc[i-1][4]*(0.000001)
        dsp_record[i][2] = vel_record[i-1][2]*(0.001)+(1/2)*anvNacc[i-1][5]*(0.000001)
        vel_record[i][0] = vel_record[i-1][0] + anvNacc[i-1][3]*(0.001)
        vel_record[i][1] = vel_record[i-1][1] + anvNacc[i-1][4]*(0.001)
        vel_record[i][2] = vel_record[i-1][2] + anvNacc[i-1][5]*(0.001)
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

def imuKalman(anvNacc):
    AccX_Value = anvNacc[:,3]#.reshape(len(anvNacc), 1)

    AccX_Variance = 0.0001
    dt = 0.001
    # transition_matrix  
    F = [[1, dt, 0.5*dt**2], 
     [0,  1,       dt],
     [0,  0,        1]]
    # observation_matrix 
    H = [0, 0, 1] 
    # transition_covariance 
    Q = [[0.2,    0,      0], 
     [  0,  0.1,      0],
     [  0,    0,  10e-4]]
    # observation_covariance 
    R = AccX_Variance
    # initial_state_mean
    X0 = [0,
      0,
      AccX_Value[0]]
    # initial_state_covariance
    P0 = [[  0,    0,               0], 
      [  0,    0,               0],
      [  0,    0,   AccX_Variance]]
    
    n_timesteps = len(AccX_Value)
    n_dim_state = 3
    filtered_state_means = np.zeros((n_timesteps, n_dim_state))
    filtered_state_covariances = np.zeros((n_timesteps, n_dim_state, n_dim_state))

    kf = KalmanFilter(transition_matrices = F, 
                  observation_matrices = H, 
                  transition_covariance = Q, 
                  observation_covariance = R, 
                  initial_state_mean = X0, 
                  initial_state_covariance = P0)

    for t in range(n_timesteps):
        if t == 0:
            filtered_state_means[t] = X0
            filtered_state_covariances[t] = P0
        else:
            filtered_state_means[t], filtered_state_covariances[t] = (
            kf.filter_update(
                filtered_state_means[t-1],
                filtered_state_covariances[t-1],
                AccX_Value[t]
            )
        )
    print(filtered_state_means[len(filtered_state_means)-1])

    #plt.ylim(-0.2, 0.2)
    #plt.ylim(-9, -10.6)
    plt.figure()
    plt.title('(Kalman)Displacement')
    plt.plot(newrec,color='red',label='x_rec')
    plt.legend()
    plt.xlabel('(ms)')
    plt.ylabel('(m)')
    plt.show()

    #f, axarr = plt.subplots(3, sharex=True)
    '''
    #axarr[0].plot(AccX_Value, label="Input AccX")
    axarr[0].plot(filtered_state_means[:, 0], "r-", label="Estimated AccX")
    axarr[0].set_title('Acceleration X')
    #axarr[0].grid()
    axarr[0].legend()
    #axarr[0].set_ylim([-4, 4])
    
    axarr[1].plot(RefVelX, label="Reference VelX")
    axarr[1].plot(filtered_state_means[:, 1], "r-", label="Estimated VelX")
    axarr[1].set_title('Velocity X')
    axarr[1].grid()
    axarr[1].legend()
    #axarr[1].set_ylim([-1, 20])

    axarr[2].plot(RefPosX, label="Reference PosX")
    axarr[2].plot(filtered_state_means[:, 0], "r-", label="Estimated PosX")
    axarr[2].set_title('Position X')
    axarr[2].grid()
    axarr[2].legend()
    #axarr[2].set_ylim([-10, 1000])
    plt.show()
    '''
def tryKalman(acc):
    dt = 0.001
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

    anvNacc = read() #data input
    aver_sample_num = 30000 #1000 sample as 1 sec.
    anvNacc = correctCoordinate(anvNacc, aver_sample_num)

    ########## Filter Field ##########
    tr=Trace()
    tr.stats['sampling_rate']=1000
    for i in range(6):
        tr.data=anvNacc[:,i]
        tr.filter("bandpass", freqmin =0.1,freqmax = 100)
        #tr.filter("highpass", freq =0.01)
        #tr.filter("lowpass", freq=100)
        anvNacc[:,i] = tr.data
    
    ########## Block of Computing Velocity & Displacement ##########

    record = tryKalman(anvNacc[:,3])
    showPlot_1_Axes(record[:,0], "Displacement", "m")

if __name__ == '__main__':
    main()