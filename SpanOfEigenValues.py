# Bentley was here
# Bentley was here again

from numpy import zeros, linalg, linspace, array, shape
from sys import exit # Terminates the code, used if n is odd
import matplotlib.pyplot as plt
from time import time # Allows run time to be found
start_time = time() # records time at start

n = 10 #number of sites
if n%2 != 0:
    print('Error: number of sites is not even, no meaningful half filling'); exit()

MinValue = 0
MaxValue = 1
NumOfPoints = 10 # Number of sample points
variable = []

'''Choose a variable - only one!!! - leave the other two commented out'''

'''Vary Ɛ'''
variable.append('epsilon')
g = 1
M = 0
epsilonArray = linspace(MinValue, MaxValue, NumOfPoints)

'''Vary M'''
# variable.append('M')
# g = 1
# MArray = linspace(MinValue, MaxValue, NumOfPoints)
# epsilon = 1

'''Vary g'''
# variable.append('g')
# gArray = linspace(MinValue, MaxValue, NumOfPoints)
# M = 0
# epsilon = 1

eigValues = []
allowedStates = [] #list of allowed states
stateLib = {} #tells you what the binary representation of each state is
stateQlib = {} #tells you what the Q at each site is, for each state
stateElib = {} #tells you what the E at each link is, for each state
pairs = [] # Will be a list containing the connected pair states

for i in range(2**(int(n/2))-1, 2**n - int(2**(n/2)) +1): #Generates binary numbers from 2^(n/2)-1 up to 2^n - 2^(n/2), (min and max half filling) all of length n
    binaryNumber = i, [int(k) for k in "{0:0{1}b}".format(i, n)] # Creates tuple, the first element is the number of the state, the second is the binary representation in list form
    if sum(binaryNumber[1]) == n/2: #Adds all 1/2 filling states to the list of allowedStates, and defines what they are in stateLib
        allowedStates.append(binaryNumber[0])
        newState = {binaryNumber[0]:binaryNumber[1]}
        stateLib.update(newState)

N = len(allowedStates) #number of states
H = zeros(shape=(N,N)) # Hamiltonian matrix, for now entirely made up of 0s
def Qfinder(state): #generates a list of the Q values at each site, for a given state, then attatches that list to the state in stateQlib
    Qlist = []
    for i in range(n):
        if i%2 == 0:
            Qlist.append(stateLib.get(state)[i])
        else:
            Qlist.append(stateLib.get(state)[i]-1)
    stateQlib.update({state:Qlist})

def Efinder(state): #generates a list of E values at each link, for a given state, then attatches that list to the state in stateElib
    Elist = []
    for i in range(n-1):
        newE = stateQlib.get(state)[i]
        if i != 0:
            newE += Elist[i-1]
        Elist.append(newE)
    stateElib.update({state:Elist})

def DiagFinder(M,g): # Generates diagonal elements of H
    for i,state in zip(range(N), allowedStates):
        Hf = M * sum([(-1) ** i * stateLib.get(state)[i] for i in range(len(stateLib.get(state)))])
        Hb = g**2 * 1/2 * sum([stateElib.get(state)[i]**2 for i in range(len(stateElib.get(state)))]) #1/2 g * the sum of E^2 on each link for a given state
        H[i][i] = Hb + Hf

def CoupleFinder(): # Finds coupled states via left motion of fermion, adds them to the list 'pairs'
    for state in allowedStates:
        for a in range(2**n - 2*int(2**(n/2)) +1):
            if state + 2**a in allowedStates:
                pairs.append([state,state+2**a])

def OffDiagFinder(epsilon): # Finds off diagonal elemenst
    for p in pairs: # Goes through the list pairs, finds the number of the paired states, and sets that element of H to e
        H[allowedStates.index(p[0])][allowedStates.index(p[1])] = epsilon
        H[allowedStates.index(p[1])][allowedStates.index(p[0])] = epsilon

def Start1(): # Runs the program for variable epsilon
    for state in allowedStates: # Runs Qfinder and Efinder for all states
        Qfinder(state)
        Efinder(state)

    DiagFinder(M,g)
    CoupleFinder()
    for epsilon in epsilonArray: # Runs OffDiagFinder for diferent values of epsilon and finds the eigenvalues of each H that creates
        OffDiagFinder(epsilon)
        eigValue = linalg.eig(H)[0] # The eigenvalues of H
        # eigVectors  = linalg.eig(H)[1] # The eigenvectors of H
        eigValues.append(sorted(eigValue)) # List of eigenvalues of H, 1 sublist per value of epsilon

    for i in range(len(eigValues[0])):
        plt.plot(epsilonArray, array(eigValues)[:,i], label=f'λ{i+1}')
    plt.title(f'Graph of EigenValues against Ɛ, for half filling of {n} sites, with M = {M} and g = {g}, using {NumOfPoints} sample points')
    plt.xlabel('Value of Ɛ')
    plt.ylabel('EigenValue')
    # plt.legend()
    # plt.xlim(0.0233333,0.0233334)
    plt.show()

def Start2(): # Runs the program for variable M
    for state in allowedStates: # Runs Qfinder and Efinder for all states
        Qfinder(state)
        Efinder(state)

    CoupleFinder()
    OffDiagFinder(epsilon)
    for M in MArray:
        DiagFinder(M,g)
        eigValue = linalg.eig(H)[0] # The eigenvalues of H
        # eigVectors  = linalg.eig(H)[1] # The eigenvectors of H
        eigValues.append(sorted(eigValue)) # List of eigenvalues of H, 1 sublist per value of epsilon

    for i in range(len(eigValues[0])):
        plt.plot(MArray, array(eigValues)[:,i], label=f'λ{i+1}')
    plt.title(f'Graph of EigenValues against M, for half filling of {n} sites, with Ɛ = {epsilon} and g = {g}, using {NumOfPoints} sample points')
    plt.xlabel('Value of M')
    plt.ylabel('EigenValue')
    # plt.legend()
    # plt.xlim(0.0233333,0.0233334)
    plt.show()

def Start3(): # Runs the program for variable g
    for state in allowedStates: # Runs Qfinder and Efinder for all states
        Qfinder(state)
        Efinder(state)

    CoupleFinder()
    OffDiagFinder(epsilon)
    for g in gArray:
        DiagFinder(M,g)
        eigValue = linalg.eig(H)[0] # The eigenvalues of H
        # eigVectors  = linalg.eig(H)[1] # The eigenvectors of H
        eigValues.append(sorted(eigValue)) # List of eigenvalues of H, 1 sublist per value of epsilon

    for i in range(len(eigValues[0])):
        plt.plot(gArray, array(eigValues)[:,i], label=f'λ{i+1}')
    plt.title(f'Graph of EigenValues against g, for half filling of {n} sites, with Ɛ = {epsilon} and M = {M}, using {NumOfPoints} sample points')
    plt.xlabel('Value of g')
    plt.ylabel('EigenValue')
    # plt.legend()
    # plt.xlim(0.0233333,0.0233334)
    plt.show()


if len(variable) == 1:
    if variable == ['epsilon']:
        Start1()
    elif variable == ['M']:
        Start2()
    elif variable == ['g']:
        Start3()

elif len(variable) >1:
    print('Error: Too many variables, comment out all but 1, see line 15')
else:
    print('Error: No variables set, comment out all but 1, see line 15')
