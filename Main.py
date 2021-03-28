# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 20:49:25 2021
@author: dbric
"""

"""
Declan Brick
Matthew Bethea
2/17/21
MAE 343
Project 1
"""

#import here
import math
#import matplotlib as m


"""
This is the first Function or option
(3 inputs | 1 output)
------------------------------------------------------------------------
"""
def mach(gamma, Beta, Theta):
    
    B = math.radians(Beta)
    T = math.radians(Theta)
    g = gamma
    
    RMach = math.sqrt((-2*math.tan(T)-2*(1/(math.tan(B)))/(g*math.tan(T)+math.tan(T)*math.cos(2*B)-2*math.sin(B)*math.sin(B)*(1/(math.tan(B))))))
    
    return RMach


"""
This is the second Function or option
(3 inputs | 1 output)
------------------------------------------------------------------------
"""
def theta(gamma, Mach, Beta):
    
    g = gamma
    M = Mach
    B = Beta
    
    """
    Calculates flow deflection angle using Equation 1 in the Project Description Parameters
    ----------
    g : Float
        Ratio of specific heats
    M : Float
        Free Stream Mach Number before Shock
    Beta : Float
        Oblique Shock Angle, in degrees
        
    Returns
    -------
    Theta: Float
        Flow Deflection Angle, in degrees
    """

    
    #Just writing rhs of Equation 1
    Tan_Theta=2*(1/math.tan(math.radians(B)))*(M**2*(math.sin(math.radians(B)))**2-1)/(M**2*(g+math.cos(math.radians(2*B)))+2)
    
    #Now grabbing the actual angle
    RTheta=math.degrees(math.atan(Tan_Theta))
    
    return RTheta
    #this was verified using a Ti-84 to output the correct answer, but many want to change output decimal precision


"""
This is the third Function or option
(3 inputs | 1 output)
------------------------------------------------------------------------
"""
def beta(gamma, Mach, Theta):
    
    
    g = gamma
    M = Mach
    T = Theta
    
    #using the functions for finding Beta max (the maximum shock angle) supported by the maximum Theta (the maximum flow deflection)
    f = lambda B: math.atan(2*(1/math.tan(math.radians(B)))*(M**2*(math.sin(math.radians(B)))**2-1)/(M**2*(g+math.cos(math.radians(2*B)))+2))-math.radians(T)
    
    g = gamma
    M = Mach
    T = Theta
    
    RBeta = []
    
    BM = M_Beta(g,M)
    
    Theta_max = M_Theta(g, M, BM)
    
    Mach_angle = Mach_a(M)

        #Strong Root for Beta function
        #This root will be between 90 and the theta max
        #F = the function
        #sA = Theta max
        #sB = 90
        #iteration is set automatically to 100

    #setting bounds
    for c in range(1,3):
        if c == 1:
            Old_A = math.degrees(BM) #p0
            Old_B = float(90) #p1
            
        if c == 2:
            Old_A = Mach_angle
            Old_B = math.degrees(BM)
        
        #verifying
##        if f(Old_A)*f(Old_B) >= 0:
##            print("failed verification")
##            return None

        #iterating to find the solution for beta
        for n in range(1, 100):
            
            new = Old_B - f(Old_B)*(Old_B - Old_A)/(f(Old_B)-f(Old_A))
            solvenew = f(new)
            
            if abs(new-Old_B) <0.0001:
                sBeta = Old_A - f(Old_A)*(Old_B - Old_A)/(f(Old_B) - f(Old_A))
                
                RBeta.append([sBeta,n])
                break

            elif n==99:
                print('failed')
            else:
                Old_A=Old_B
                Old_B=new
   
    return RBeta


"""
This is a sub function of function 3
(2 inputs | 1 output)
---------------------------
"""
def M_Beta(g,M):
    #print("veriables")
    #print(g)
    #print(M)
    f1 = 1+((g-1)/2)*(M**2)+((g+1)/16)*(M**4)
    f2 = math.sqrt((g+1)*f1)
    f3 = (((g+1)/4)*(M**2)+f2-1)
    f4 = math.sqrt((1/(g*(M**2)))*f3)
    BM = math.asin(f4)
    #print("Beta Max")
    #print(BM)
    return BM


"""
This is a sub function of function 3
(3 inputs | 1 output)
---------------------------
"""
def M_Theta(gamma, Mach, Beta):

    #print("Varribles")
    #print(gamma)
    #print(Mach)
    #print(Beta)
    
    g = gamma
    M = Mach
    B = Beta
    
    top = ((M**2)*(math.sin(B)*math.sin(B))-1)*(1/math.tan(B))
    bot = (((1/2)*(g+1)*(M**2))-((M**2)*(math.sin(B)*math.sin(B)))+1)
    
    TM = math.degrees(1/(math.tan(top/bot)))

    #print("The maximum flow deflection is ")
    #print(TM)
    
    return TM


"""
This is a sub function of function 3
(1 inputs | 1 output)
---------------------------
"""
def Mach_a(Mach):
    
    M = Mach
    
    mu = math.degrees(math.asin(1/M))

    #print(mu)
    
    return mu

"""
This section is about graphing
--------------------------------

def Cp(M,g):
    for i in range(0,len(M)):
        plot=[]
        m1=M[i]
        CP=[]
        Beta=[]
        Theta=[]
        for ii in range(0,21): #solve for Beta values
            Theta.append(ii)
            Beta.append(beta(g,m1,ii)[1][0])
    
        for ii in range(0,21): #solve for pressure values
            print(Beta[ii])
            Mn1=m1*math.sin(math.radians(Beta[ii]))
            Pressure_Ratio=1+(Mn1**2-1)*(2*g/(g+1))
            CP.append(2/(g*m1**2)*(Pressure_Ratio-1))
        print('Line'+str(i))
        m.pyplot.plot(Theta,CP, label='M='+str(m1))
    m.pyplot.title('Pressure Coefficient vs. Flow Deflection Angle')
    m.pyplot.xlabel('Flow Deflection Angle (°)')
    m.pyplot.ylabel('Pressure Coefficient')
    m.pyplot.figlegend(loc='best')
    return 
"""



"""
                        Main Code area
                        This is the start of the main code
"""

"""
                        This matrix is the options or function dictionary
                        Left side are the keys that gets printed to the user
                        The right is the fuction to be called when the user matches with
"""
functionMatrix = {
  "Mach":mach,
  "Theta":theta,
  "Beta":beta ,
}

"""
                        This matrix is a dictionary of the  varriables required in each function choice
                        The left side is the keys that match with the fuction
                        The right side are arrays/matrix that hold the varribles (order matters | list = fuction order)
"""
lst = {
    "Mach":["gamma", "Beta(°)", "Theta(°)"],
       "Theta":["gamma", "Mach", "Beta(°)"],
       "Beta":["gamma", "Mach", "theta(°)"]
    }

ver = []
results = []

"""
                        This print fuction prompts the user with the avalible functions
"""
print("Available Functions:")           # Header
for function in functionMatrix.keys():  # basically a for loop that uses creates varrible called functions that it stores the current key in
  print(function)                       # Prints the keys while the for loop runs

"""
                        Ask the user to choose what they would like to do
"""
chosen = str(input("Which Function Would You Like to run? "))   #Prompts for the users input and stores it
print(chosen)                                                   #Shows the user whick option they choose

"""
                        Uses the users choice to prompt them with which varribles are required
"""
variableList = lst[chosen]                          #Based on the users choice it matches it with the lst dictionary key and stores the
                                                    # specific sub matrix in its own new list/matrix
for var in variableList:                            #Creates a for loop that runs for the number of varribles in the VariableList matrix
                                                    # for every instance it saves the current varrible as var that iteration of the loop
    ele = float(input(f"Please input {var}: "))     # prompts the user for the numeric value of specific variables 
    ver.append(ele)                                 # This line saves the  value the value the user entered and appends it to a 
                                                    # an empty list in the order they are asked in the function

""" This is where the magic happens"""
results = functionMatrix[chosen](*ver)      # This line runes the function chosen by the user with the variables required by the function
a = results                                        
                   
"""
Should be needed if the fuction needs to output to variables will work on it a different way
-----

if isinstance(results, list):
    a = results[0]
    b = results[1]

    # list returned
else:
    a = results
    # 1 value returend
"""

"""
                                        Here we are printing out the results to the user
                                        We started with the specific information to each seaction
                                        Then the common information is done last
"""
if chosen == "Mach":
    print("**************Results**************\n\n")
    print("----- Given -----")
    print("Shock Angle: ", ver[1])
    print("Flow Deflection angle: ", ver[2])
    print("")
    print("----- Results -----")
    print("Mach Number: ", round(a,4))

    gamma = ver[0]
    M = a
    theta = ver[2]
    Beta = ver[1]

if chosen == "Theta":
    print("**************Results**************\n\n")
    print("----- Given -----")
    print("Mach Number: ", ver[1])
    print("Shock Angle: ", ver[2])
    print("")
    print("----- Results -----")
    print("Flow Deflection angle: ", round(a,4))

    gamma = ver[0]
    M = ver[1]
    theta = a
    Beta = ver[2]
    
if chosen == "Beta":
    print("**************Results**************\n\n")
    print("----- Given -----")
    print("Mach Number: ", ver[1])
    print("Flow Deflection angle: ", ver[2])
    print("")
    print("----- Results -----")
    print("# of Iterations to Convergence", a[1])
    print("Shock Angle", round(a,4))
    
    gamma = ver[0]
    M = ver[1]
    theta = ver[2]
    Beta = a[0]
    
    
    
Mn1=M*math.sin(math.radians(Beta))
Mn2=((Mn1**2+(2/(gamma-1)))/((2*gamma/(gamma-1)*Mn1**2)-1))**0.5
M2=Mn2/(math.sin(math.radians(Beta-theta)))
Pressure_Ratio=1+(Mn1**2-1)*(2*gamma/(gamma+1))
Density_Ratio=((gamma+1)*Mn1**2)/((gamma-1)*Mn1**2+2)
Temperature_ratio=Pressure_Ratio/Density_Ratio
print("M1 normal: ", round(Mn1,4))
print("M2 normal: ", round(Mn2,4))
print("M2: ", round(M2,4))
print("P2/P1: ", round(Pressure_Ratio,4))
print("T2/T1: ",round(Temperature_ratio,4))



"""
Mn1=M*math.sin(math.radians(Beta))
Mn2=((Mn1**2+(2/(gamma-1)))/((2*gamma/(gamma-1)*Mn1**2)-1))**0.5
M2=Mn2/(math.sin(math.radians(Beta-Theta)))
Pressure_Ratio=1+(Mn1**2-1)*(2*gamma/(gamma+1))
Density_Ratio=((gamma+1)*Mn1**2)/((gamma-1)*Mn1**2+2)
Temperature_ratio=Pressure_Ratio/Density_Ratio
"""

#Matt is Awesome
