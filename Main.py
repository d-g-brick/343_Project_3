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
    

    





def TM(g, M, ThetaS, dTheta):
    """
    Solves Taylor-Maccoll Eq for Cone angle given shock angle and Mach

    Parameters
    ----------
    g : Float
        Ratio of Specific Heats
    M : Float
        Free-Stream Mach Number
    ThetaS : Float
        Shock Wave Angle, in Deg

    Returns
    -------
    None.

    """

    #First Solve for the flow deflection angle and M_2 directly behind the shock
    # Using Oblique Shock Relations 
    # deflection
    delta=theta(g,M,ThetaS) #calling function developed in Project 1
    #M_2
    Mn1=M*math.sin(math.radians(ThetaS))
    Mn2=((Mn1**2+(2/(g-1)))/((2*g/(g-1)*Mn1**2)-1))**0.5
    M2=Mn2/(math.sin(math.radians(ThetaS-delta)))
    #Now Grab V' using Mach (Eq 10.16 in book)
    Vp=1/math.sqrt(1+(2/((g-1)*M2**2)))
    #correct so far
    #Now Grab V'_r and V'_theta from geometry (Figure 10.4 in book)
    A=ThetaS-delta
    Angle=math.radians(A)
    Vpr=Vp*math.cos(Angle)
    Vpt=-Vp*math.sin(Angle) 
    #checked with Kopal's tables, and these values look roughly correct
    
    
    #Runge-Kutta time, making a subfunction for this
    [Theta,Vr,Vt,it]=TM_RK(g,Vpr,Vpt,ThetaS,dTheta) #think 0.01 is a reasonable deltaTheta
    
    return Theta,Vr,Vt,Mn1,M2,it
    



def TM_RK(g,Vpr,Vpt, ThetaS, dTheta):
    """
    

    Parameters
    ----------
    g : Float
        Ratio of Specific Heats
    Vpr : Float
        Radial velocity just behind shock
    Vpt : Float
        Velocity in Angular direction just behind shock
    ThetaS: Float
        Shock Angle
    dTheta: Float
        Angle increment, in deg, should be positive


    Returns
    -------
    None.

    """
    #note that the first equations are solving for Vr and that dr/dtheta=Vtheta
    
    maxiter=ThetaS/-dTheta #this is the number of iteratoins to step from 90 degrees (normal shock) to 0, if no solution is found, you're not going to
    it=0#initialize iteriation couinter
    theta=ThetaS #first angle
    
    rdt=math.radians(dTheta)#thinking this might be an issue
    #need a few arrays to return
    Tarr=[theta]#array of angles
    Vrarr=[Vpr]#array of radial velocities
    Vtarr=[Vpt]#array of angular velocities
    
    #writing the f_1 and f_2 equations as a subfunction below
    #angle starts at thetas and then subtracts dtheta
    
    #loop time
    
    while ((it<maxiter) & (Vpt <0.0)): #not sure about the Vpt condition, need further testing    
        
        tstep=theta+(dTheta/2)
        
        #RK k values
        k11=rdt*Vpt
        k12=rdt*f2(g,theta,Vpr,Vpt)
        k21=rdt*(Vpt+(k12/2))
        k22=rdt*f2(g,tstep,Vpr+(k11/2),Vpt+(k12/2))
        k31=rdt*(Vpt+(k22/2))
        k32=rdt*f2(g,tstep,Vpr+(k21/2),Vpt+(k22/2))
        k41=rdt*(Vpt+k32)
        k42=rdt*f2(g,tstep,Vpr+k31,Vpt+k32)
        
        #find the addition to the velocities
        changeR=(k11+2*k21+2*k31+k41)/6
        changeT=(k12+2*k22+2*k32+k42)/6
        
        #add the addition
        Vpr=Vpr+changeR
        Vpt=Vpt+changeT 
        
        #store in array
        Vrarr.append(Vpr)
        Vtarr.append(Vpt)
        
        #calculate theta
        theta=theta+dTheta #working backwards
        Tarr.append(theta)#add ther angle to the array
        
        
        it+=1#count iteration
        
    return Tarr,Vrarr,Vtarr,it

def f2(g, theta,Vr,Vt):
    """
    

    Parameters
    ----------
    g : Float
        Ratio of Specific Heats
    theta : Float
        Angle, in deg 
    Vr : Float
        Radial Velocity
    Vt : Float
        Angular velocity

    Returns
    -------
    f2, derivative of the angular velocity wrt angle.  

    """
    num=-(((g-1)/2)*(1-Vr**2-Vt**2)*(2*Vr+(Vt/math.tan(math.radians(theta)))))+(Vt**2*Vr)
    denom=(((g-1)/2)*(1-Vr**2-Vt**2))-Vt**2
    f2=num/denom
    return f2
    

def TM_Shock(g,M,ThetaC):
    """
    Solves Taylor-Maccoll Equation for shock angle given Mach and Cone angle
    

    Parameters
    ----------
    g : Float
        Ratio of Specific Heats
    M : float
        Mach number
    ThetaC : float
        cone anlge, in Deg

    Returns
    -------
    Shock angle

    """
    
    #First define the function we're trying to find the zero of
    f= lambda SG: TM(g,M,SG)[0][-1]-ThetaC #want the determined cone angle to be the real cone angle, the indexing should grab the last element of the angles, ie: the calculated cone angle
    
    
    #define the bounds
    p0=ThetaC#can't be smaller than the cone
    p1=75 #from PLot in Lecture 16, eye ball for a reasonable threshold. Note that this won't work for M<1.1 
    
    
    for n in range(1, 100):
        print(n)
        q0=f(p0)
        q1=f(p1)
        p= p1 - q1*((p1 - p0)/(q1-q0))
        
        if abs(p-p1)<1e-6:
            
            ThetaS=p
            break

        else:
            p0=p1
            p1=p
    return ThetaS
    
    
    
    
    
    
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
                        =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                        =-=-=-=-=-=  Main Code  =-=-=-=-=-=
                        =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
"""

"""This matrix is the options or function dictionary
Left side are the keys that gets printed to the user"""

chosen = 0
ver = []
results = []
cheese =[]

functionMatrix = {
  "(1) Cone Angle":TM,
  "(2) Shock Angle":TM_Shock
}

"""
                        This matrix is a dictionary of the  varriables required in each function choice
                        The left side is the keys that match with the fuction
                        The right side are arrays/matrix that hold the varribles (order matters | list = fuction order)
"""
lst = {
    "1":["gamma", "Mach", "Shock Angle", "dtheta"],
    "2":["gamma", "Mach", "Cone Angle", "dtheta"]
}


"""This print fuction prompts the user with the avalible functions"""

print("Available Functions:")           # Header
for function in functionMatrix.keys():  # basically a for loop that uses creates varrible called functions that it stores the current key in
  print(function)                       # Prints the keys while the for loop runs


"""Ask the user to choose what they would like to do"""

chosen = input("Which Function Would You Like to run? ")   #Prompts for the users input and stores it

if chosen == "1":
    func = "(1) Cone Angle"
    print("Chosen: Cone Angle")
    
if chosen == "2":
    func = "(2) Shock Angle"
    print("Chosen: Shock Angle")

"""Uses the users choice to prompt them with which varribles are required"""

variableList = lst[chosen]                          #Based on the users choice it matches it with the lst dictionary key and stores the
                                                    # specific sub matrix in its own new list/matrix
for var in variableList:                            #Creates a for loop that runs for the number of varribles in the VariableList matrix
                                                    # for every instance it saves the current varrible as var that iteration of the loop
    ele = float(input(f"Please input {var}: "))     # prompts the user for the numeric value of specific variables 
    ver.append(ele)                                 # This line saves the  value the value the user entered and appends it to a 
                                                    # an empty list in the order they are asked in the function

""" This is where the magic happens"""
results = functionMatrix[func](*ver)      # This line runes the function chosen by the user with the variables required by the function
                   

"""Here we are printing out the results to the user
We started with the specific information to each seaction
Then the common information is done last"""

print("\n --- Inital Contions --- ")
print("Gamma: ", ver[0])
print("Mach: ", ver[1])
if chosen == "1":
    print("Shock: ", ver[2])
if chosen == "2":
    print("Cone Angle: ",ver[2])
print("dtheta: ", ver[3])

if chosen == "1":
    M2 = results[4]
    Mn1 = results[3]
    
if chosen == "2":
    cheese = TM(ver[0], ver[1], results, ver[3])
    M2 = cheese[4]
    Mn1 = cheese[3]
    
Pr=1+(Mn1**2-1)*(2*ver[0]/(ver[0]+1))
Dr=((ver[0]+1)*Mn1**2)/((ver[0]-1)*Mn1**2+2)
Tr=Pr/Dr

print("\n --- Conditions behind shock --- ")
print("Mach: ", round(M2, 4))
print("Pressure: ", round(Pr, 4))
print("Temp: ", round(Tr, 4))
print("Density: ", round(Dr, 4))


print("\n --- Final Conditions ---")
if chosen == "1":
    a = results[0]
    b = results[1]
    c = results[2]
    
    print("Shock Angle: ", a[0])
    print("Cone Angle: ", round(a[-1], 4))
    print("V_r: ", round(b[-1], 4))
    print("V_Theta: ", round(c[-1], 4))
    
    V = math.sqrt(((b[-1])**2) + ((c[-1])**2))
    
    print("V': ", round(V, 4))

    print("Iterations to answer: ", results[5])