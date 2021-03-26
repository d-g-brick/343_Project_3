# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 14:00:00 2021
@author: dbric
"""

"""
Declan Brick
Matthew Bethea
2/17/21
MAE 343
Project 2
"""

#import here
import math
import matplotlib as m

# Project 1 Code needed for Shock, only need to calculate Beta as given M_1 and Theta
# also need auxillarly functions from Project 1

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


def Mach_a(Mach):
    
    M = Mach
    
    mu = math.degrees(math.asin(1/M))

    #print(mu)
    
    return mu

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


    #setting bounds
    for c in range(1,3):
        if c == 1:
            Old_A = math.degrees(BM) #p0
            Old_B = float(90) #p1
            
        if c == 2:
            Old_A = Mach_angle
            Old_B = math.degrees(BM)

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


# Project 2 code 

def Nu(gamma, Mach):
    
    M = Mach
    M2 = M**2
    g = gamma
    gp1 = g+1
    gm1 = g-1
    
    VofM = math.degrees(math.sqrt(gp1/gm1)*math.atan(math.sqrt((gm1/gp1)*(M2-1)))-math.atan(math.sqrt(M2-1)))
    
    return VofM



def M2(gamma, Nu):
    '''

    Parameters
    ----------
    gamma : Ratio of specific heats, assumed constant.
    Nu : Prandtl Meyer value for M_2, expect input in degrees

    Returns
    -------
    M_2 : Mach number 

    '''    
    #shorthands
    g = gamma
    N= Nu
    
    #the Prandtl Meyer function-the value, should be ~0 when M_2 is given
    f = lambda M: math.sqrt((g+1)/(g-1))* math.atan(math.sqrt((g-1)*(M**2-1)/(g+1)))-math.atan(math.sqrt(M**2-1))-math.radians(N)
    



    #setting bounds
    # limits choosen based on limits of Table A.5
    Old_A = 1 #p0
    Old_B = 7.3 #p1
            
        #iterating to find M
    for n in range(1, 100):
            '''
            Old Troubleshooting Code
            print(n)
            print(Old_A)
            print(Old_B)
            print('-----------------------')
            '''
            q_A=f(Old_A)#value of f at first bound
            q_B=f(Old_B)#value of f at second bound
            new = Old_B - q_B*(Old_B - Old_A)/(q_B-q_A)#secant equation for finding new guess
            
            if abs(new-Old_B) <0.00001:#if within convergence tolerance
                M_2= Old_A - f(Old_A)*(Old_B - Old_A)/(f(Old_B) - f(Old_A))
                break

            elif n==99:#wayyyyy too many iterations
                print('failed')
            else:
                #set new bounds
                Old_A=Old_B
                Old_B=new
                #get imaginary number issues if guess drops below M=1
                if Old_A<1:
                    Old_A=1
                if Old_B<1:
                    Old_B=1
   
    return M_2




#Main Code area


#this should solve for the top half of an airfoil, assumes straight line sides


ver=[]

variableList = ["gamma", "Mach", "alpha(°)","Epsilon(°)", "Theta(°)"]
for var in variableList:
    ele = float(input(f"Please input {var}: "))
    ver.append(ele)

#do solving here
if ver[2]<ver[3]: #have shock then expansion
    gamma=ver[0]
    theta=ver[3]-ver[2]# deflection is epsilon-alpha
    Beta=beta(gamma,ver[1],theta)
    Beta=Beta[1][0]
    # determining ratios and values
    Mn1=ver[1]*math.sin(math.radians(Beta))
    Mn2=((Mn1**2+(2/(gamma-1)))/((2*gamma/(gamma-1)*Mn1**2)-1))**0.5
    M_2=Mn2/(math.sin(math.radians(Beta-theta)))
    Pressure_Ratio_1=1+(Mn1**2-1)*(2*gamma/(gamma+1))
    Density_Ratio_1=((gamma+1)*Mn1**2)/((gamma-1)*Mn1**2+2)
    Temperature_ratio_1=Pressure_Ratio_1/Density_Ratio_1
    
    # now for expansion
    Nu_M2=Nu(gamma,M_2)
    Mu2=math.degrees(math.asin(1/M_2))
    Nu_M3=Nu_M2+ver[4]
    M_3=M2(gamma, Nu_M3)
    Mu3=math.degrees(math.asin(1/M_3))
    
    Temperature_ratio_2=(1+(gamma-1)*M_3**2/2)
    Pressure_Ratio_2=Temperature_ratio_2**(gamma/(gamma-1))
    
    
    
    print("**************Results**************\n\n")
    print("----- Given -----")
    
    print("Ratio of Specific Heats: ", ver[0])
    print("Initial Mach: ",ver[1])
    print("Angle of Attack: ", ver[2])
    print("Half Angle: ", ver[3])
    print("Expansion Deflection Angle:", ver[4])
    print("Based on angle of attack and half angle, Region 1 is upstream,")
    print("there is an oblique shock between Region 1 and Region 2, and  ")
    print("there is an Prandtl-Meyer expansion wave between Region 2 and Region 3.")
    print("----- Results -----")
    print("Flow Deflection Angle for Shock Wave: ", theta,"(°)")
    print("P2/P1: ", round(Pressure_Ratio_1,4))
    print("T2/T1: ",round(Temperature_ratio_1,4))
    print("Mach Number in Region 2: ", round(M_2,4))
    print("Nu2: ", round(Nu_M2,4),"(°)")
    print("mu2: ", round(Mu2,4),"(°)")
    print("Mach Number in Region 3: ", round(M_3,4))
    print("Nu3: ", round(Nu_M3,4),"(°)")
    print("mu3: ", round(Mu3,4),"(°)")
    print("P3/P2: ", round(Pressure_Ratio_2,4))
    print("T3/T2: ",round(Temperature_ratio_2,4))
    print("P3/P1: ", round(Pressure_Ratio_2*Pressure_Ratio_1,4))
    print("T1/T1: ",round(Temperature_ratio_2*Temperature_ratio_1,4))
    
    
    
elif ver[2]==ver[3]:#just have expansion
    gamma=ver[0]
    
    # now for expansion
    Nu_M1=Nu(gamma,ver[1])
    Mu1=math.degrees(math.asin(1/ver[1]))
    Nu_M2=Nu_M1+ver[4]
    M_2=M2(gamma, Nu_M2)
    Mu2=math.degrees(math.asin(1/M_2))
    
    Temperature_ratio=(1+(gamma-1)*M_2**2/2)
    Pressure_Ratio=Temperature_ratio**(gamma/(gamma-1))
    
    print("**************Results**************\n\n")
    print("----- Given -----")
    
    print("Ratio of Specific Heats: ", ver[0])
    print("Initial Mach: ",ver[1])
    print("Angle of Attack: ", ver[2])
    print("Half Angle: ", ver[3])
    print("Expansion Deflection Angle:", ver[4])
    print("Based on angle of attack and half angle, Region 1 is upstream,")
    print("Region 2 has the upstream conditions, and there is an ")
    print("Prandtl-Meyer expansion wave between Region 2 and Region 3.")
    print("")
    print("----- Results -----")
    print("Nu1: ", round(Nu_M1,4))
    print("mu1: ", round(Mu1,4))
    print("Mach Number in Region 2: ", round(M_2,4))
    print("Nu2: ", round(Nu_M2,4))
    print("mu2: ", round(Mu2,4))
    print("P2/P1: ", round(Pressure_Ratio,4))
    print("T2/T1: ",round(Temperature_ratio,4))
    
elif ver[2]>ver[3]:#have expansion then expansion
    gamma=ver[0]
    
    # now for expansion 1
    Nu_M1=Nu(gamma,ver[1])
    Mu1=math.degrees(math.asin(1/ver[1]))
    Nu_M2=Nu_M2+(ver[2]-ver[3])
    M_2=M2(gamma, Nu_M3)
    Mu2=math.degrees(math.asin(1/M_2))
    Temperature_ratio_1=(1+(gamma-1)*M_2**2/2)
    Pressure_Ratio_1=Temperature_ratio_2**(gamma/(gamma-1))
    
    #and expansion 2
    #Nu_M2 is the same as is M_2
    Nu_M3=Nu_M2+ver[4]
    M_3=M2(gamma, Nu_M3)
    Mu3=math.degrees(math.asin(1/M_3))
    Temperature_ratio_2=(1+(gamma-1)*M_3**2/2)
    Pressure_Ratio_2=Temperature_ratio_2**(gamma/(gamma-1))
    
    print("**************Results**************\n\n")
    print("----- Given -----")
    
    print("Ratio of Specific Heats: ", ver[0])
    print("Initial Mach: ",ver[1])
    print("Angle of Attack: ", ver[2])
    print("Half Angle: ", ver[3])
    print("Expansion Deflection Angle:", ver[4])
    print("Based on angle of attack and half angle, Region 1 is upstream,")
    print("there is an Prandtl-Meyer expansion wavw between Region 1 and Region 2, and  ")
    print("there is an Prandtl-Meyer expansion wave between Region 2 and Region 3.")
    print("----- Results -----")
    print("Flow Deflection Angle for First Expansion Wave: ", ver[2]-ver[3],"(°)")
    print("P2/P1: ", round(Pressure_Ratio_1,4))
    print("T2/T1: ",round(Temperature_ratio_1,4))
    print("Mach Number in Region 2: ", round(M_2,4))
    print("Nu2: ", round(Nu_M2,4),"(°)")
    print("mu2: ", round(Mu2,4),"(°)")
    print("Mach Number in Region 3: ", round(M_3,4))
    print("Nu3: ", round(Nu_M3,4),"(°)")
    print("mu3: ", round(Mu3,4),"(°)")
    print("P3/P2: ", round(Pressure_Ratio_2,4))
    print("T3/T2: ",round(Temperature_ratio_2,4))
    print("P3/P1: ", round(Pressure_Ratio_2*Pressure_Ratio_1,4))
    print("T1/T1: ",round(Temperature_ratio_2*Temperature_ratio_1,4))
    
        
else:
    print('Error: Check Angle of attack and half angle inputs.')    



