"""Anas Lasri Doukkali, CID: 01209387
M3C 2018 Homework 3
Contains five functions:
    plot_S: plots S matrix -- use if you like
    simulate2: Simulate tribal competition over m trials. Return: all s matrices at final time
        and fc at nt+1 times averaged across the m trials.
    performance: To be completed -- analyze and assess performance of python, fortran, and fortran+openmp simulation codes
    analyze: To be completed -- analyze influence of model parameter, g
    visualize: To be completed -- generate animation illustrating "interesting" tribal dynamics
"""
import numpy as np
import time
import matplotlib.pyplot as plt
from m1 import tribes as tr #assumes that hw3_dev.f90 has been compiled with: f2py --f90flags='-fopenmp' -c hw3_dev.f90 -m m1 -lgomp
#May also use scipy and time modules as needed


def plot_S(S):
    """Simple function to create plot from input S matrix
    """
    ind_s0 = np.where(S==0) #C locations
    ind_s1 = np.where(S==1) #M locations
    plt.plot(ind_s0[1],ind_s0[0],'rs')
    plt.plot(ind_s1[1],ind_s1[0],'bs')
    plt.show()
    return None
#------------------


def simulate2(N,Nt,b,e,g,m):
    """Simulate m trials of C vs. M competition on N x N grid over
    Nt generations. b, e, and g are model parameters
    to be used in fitness calculations.
    Output: S: Status of each gridpoint at end of simulation, 0=M, 1=C
            fc_ave: fraction of villages which are C at all Nt+1 times
                    averaged over the m trials
    """
    #Set initial condition
    S  = np.ones((N,N,m),dtype=int) #Status of each gridpoint: 0=M, 1=C
    j = int((N-1)/2)
    S[j,j,:] = 0
    N2inv = 1./(N*N)

    fc_ave = np.zeros(Nt+1) #Fraction of points which are C
    fc_ave[0] = S.sum()

    #Initialize matrices
    NB = np.zeros((N,N,m),dtype=int) #Number of neighbors for each point
    NC = np.zeros((N,N,m),dtype=int) #Number of neighbors who are Cs
    S2 = np.zeros((N+2,N+2,m),dtype=int) #S + border of zeros
    F = np.zeros((N,N,m)) #Fitness matrix
    F2 = np.zeros((N+2,N+2,m)) #Fitness matrix + border of zeros
    A = np.ones((N,N,m)) #Fitness parameters, each of N^2 elements is 1 or b
    P = np.zeros((N,N,m)) #Probability matrix
    Pden = np.zeros((N,N,m))
    #---------------------

    #Calculate number of neighbors for each point
    NB[:,:,:] = 8
    NB[0,1:-1,:],NB[-1,1:-1,:],NB[1:-1,0,:],NB[1:-1,-1,:] = 5,5,5,5
    NB[0,0,:],NB[-1,-1,:],NB[0,-1,:],NB[-1,0,:] = 3,3,3,3
    NBinv = 1.0/NB
    #-------------

    #----Time marching-----
    for t in range(Nt):
        R = np.random.rand(N,N,m) #Random numbers used to update S every time step

        #Set up coefficients for fitness calculation
        A = np.ones((N,N,m))
        ind0 = np.where(S==0)
        A[ind0] = b

        #Add boundary of zeros to S
        S2[1:-1,1:-1,:] = S

        #Count number of C neighbors for each point
        NC = S2[:-2,:-2,:]+S2[:-2,1:-1,:]+S2[:-2,2:,:]+S2[1:-1,:-2,:] + S2[1:-1,2:,:] + S2[2:,:-2,:] + S2[2:,1:-1,:] + S2[2:,2:,:]

        #Calculate fitness matrix, F----
        F = NC*A
        F[ind0] = F[ind0] + (NB[ind0]-NC[ind0])*e
        F = F*NBinv
        #-----------

        #Calculate probability matrix, P-----
        F2[1:-1,1:-1,:]=F
        F2S2 = F2*S2
        #Total fitness of cooperators in community
        P = F2S2[:-2,:-2,:]+F2S2[:-2,1:-1,:]+F2S2[:-2,2:,:]+F2S2[1:-1,:-2,:] + F2S2[1:-1,1:-1,:] + F2S2[1:-1,2:,:] + F2S2[2:,:-2,:] + F2S2[2:,1:-1,:] + F2S2[2:,2:,:]

        #Total fitness of all members of community
        Pden = F2[:-2,:-2,:]+F2[:-2,1:-1,:]+F2[:-2,2:,:]+F2[1:-1,:-2,:] + F2[1:-1,1:-1,:] + F2[1:-1,2:,:] + F2[2:,:-2,:] + F2[2:,1:-1,:] + F2[2:,2:,:]

        P = (P/Pden)*g + 0.5*(1.0-g) #probability matrix
        #---------

        #Set new affiliations based on probability matrix and random numbers stored in R
        S[:,:,:] = 0
        S[R<=P] = 1

        fc_ave[t+1] = S.sum()
        #----Finish time marching-----

    fc_ave = fc_ave*N2inv/m

    return S,fc_ave
#------------------


def performance(n,nt,b,e,g,numthreads):
    """I would like to start by first saying that having a computer with 4
    physical cores I took the liberty of running all my calculations on a maximum
    of 4 threads.
    I first start by only using the simulate2_omp function and analyzing the
    speed up for different number of threads with different number of simulations
    I quickly notice the following; for low values of M(~100-200) the calculation
    speed is inversely proportional(in a linear maner) to the number of threads.
    This is the expected result as for example for the same calculation one thread
    will take 0.40sec while the 4-threads version will take 0.10sec in total.
    I would also like to point out that my function can also output the wall time
    of all the calculations for different number of threads, however I will comment
    to make the marking process easier, I found this information very useful while
    looking for trends as it led me to see the following two important things:
    1) You will notice in my first plot where I compare the speed against the number
    of threads a small jump(M~1000) in the speed of the code. This can be explained following
    what was said in the lecture about computer architecture. This jump is due to the
    amount of data being used for the calculations, as was explained in said lecture
    higher data storage units are located further away in the computer than the lower
    data storage units, from the calculations and arithmetic units.
    2) I will now like to take some time to mention Ahmdal's law from lectures,
    While looking at the parallelized code that I gave from Fortran you can notice
    that I have a section at the top that amounts to approximately 10% of the code
    which is not parallelized, This explains using the T(N) = s + p/N equation
    why as 'M' gets larger four threads run no more at four times the speed of one
    thread. This was really interesting to notice in the code and analysis of
    performance.
    I know plot one other plot with the code running for the same parameters with
    the three types of implementations. Python, as exprected, is the slowest code
    with orders of magnitude slower speed than the parallelized version and also
    as expected slower than the Fortran f2py version since we are using a compiled
    language. It is also interesting to note that this difference in speed will
    a lot of difference as 'M' gets bigger and for real life simulations of our
    code.
    """
    tr.tr_g = g
    tr.tr_b = b
    tr.tr_e = e

    #We will start by analyzing the performance of simulate_omp by itself
    m = np.linspace(0,1500,60)
    m[0] = 1
    l = len(m)
    wall_time = np.zeros((numthreads,l),dtype=float)
    wall_time1 = np.zeros((l),dtype=float)
    wall_time2 = np.zeros((l),dtype=float)
    su_time = np.zeros((numthreads,l),dtype=float)

    for j in range(l):
        for i in range(numthreads):
            tr.numthreads = i+1
            wall_time0 = time.time()
            tr.simulate2_omp(n,nt,m[j])
            wall_time[i,j] = time.time()-wall_time0
    su_time[:,:] = wall_time[0,:]/wall_time[:,:]

    for j in range(l):
        print('M:  ',m[j])
        for i in range(numthreads):
            print('Number of threads', i+1,'wall_time is :',wall_time[i,j])
            print('Number of threads', i+1,'speed up time is :',su_time[i,j])
        print()

    for j in range(l):
        wall_time0 = time.time()
        simulate2(n,nt,b,e,g,int(m[j]))
        wall_time2[j] = time.time()-wall_time0

    for j in range(l):
        wall_time0 = time.time()
        tr.simulate2_f90(n,nt,m[j])
        wall_time1[j] = time.time()-wall_time0



    for i in range(numthreads):
        plt.hold(True)
        #plt.plot(m,wall_time1)
        #plt.plot(m,wall_time2)
        plt.plot(m,wall_time[i,:])
    plt.title('Anas Lasri, CID:01209387 \n number of simuations against speed')
    plt.xlabel('Number of simulations: M')
    plt.ylabel('wall time')
    plt.legend(('one thread', 'two threads', 'three threads', 'four threads'),loc='best')
    plt.savefig('hw31.png', dpi=400)
    plt.show()

    plt.plot(m,wall_time2)
    plt.plot(m,wall_time1)
    plt.plot(m,wall_time[3,:])
    plt.title('Anas Lasri, CID:01209387 \n number of simuations againt speed')
    plt.xlabel('Number of simulations: M')
    plt.ylabel('wall time')
    plt.legend(('Python', 'f2py Fortran', 'Fortran+OpenM'),loc='best')
    plt.savefig('hw32.png', dpi=400)





    return None  #Modify as needed

def analyze(n,nt,m,e):
    """Analyze influence of model parameter, g.
    Modify the contents of the tuple, input, as needed
    When display is True, figures equivalent to those
    you are submitting should be displayed
    """
    tr.tr_e = e
    #What I will do is have 3 figures where I will plot one augmenting b and g
    #in the 2nd I will augment b and decrease g and in the last I will decrease
    #both
    lin_b = np.linspace(1.1,1.5,nt+1)
    lin_g = np.linspace(0.80,1.0,nt+1)
    fc_ave = np.zeros((len(lin_b),len(lin_g)))
    for tr.tr_b in enumerate(lin_b):
        for tr.tr_g in enumerate(lin_g):
            _,fc_ave[:,:] = tr.simulate2_omp(n,nt,m)
    for i in range(nt):
        plt.hold(True)
        plt.plot(np.linspace(1,nt,nt+1),fc_ave[:,i])


    return None #Modify as needed



def visualize():
    """Generate an animation illustrating the evolution of
        villages during C vs M competition
    """

    return None #Modify as needed


if __name__ == '__main__':
    #Modify the code here so that it calls performance analyze and
    # generates the figures that you are submitting with your code

    input_p = None
    output_p = performance(input_p) #modify as needed

    input_a = None
    output_a = performance(input_a)
