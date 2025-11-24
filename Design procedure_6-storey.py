import math
import numpy as np
import Predict_gamma
###预设
##周期
T=0.0724*(26.4)**0.8*1.4# Empirical period formula (s)
T=1.2#Input custom period (s)
print('T=',T)
M=3778571.428571428##Mass kg
#Prdefine hysteretic parameters
theta1_ = 0.3###Prdefine target drift ratio (%) (θ1,θ2,θ3,θ4)
theta2_ = 1.5
theta3_ = 1.8
theta4_ = 2.25
a1_ = 0.05#stiffness ratios α1 α2 α3
a2_ = 0.6
a3_ = 0.1
ξ_=0.05#ductility ζ
p = 0.5#Input probability performance level pne
ξ1_ = theta2_/theta1_#displacement ratios ζ1 ζ2
ξ2_ = theta3_/theta1_
miu = theta4_/theta1_
Ss = 1.0#Input the mapped acceleration parameters according to ASCE 7-22
S1 = 0.85
Fv = 1.3
Fa = 2.4
Sms = Ss*Fa
Sm1 = S1*Fv
T0 = 0.2*(Sm1/Sms)
Ts = (Sm1/Sms)
gama_list = Predict_gamma.predict(a1_,a2_,a3_,ξ_,ξ1_,ξ2_,miu,T,p)
gama_DBE  = gama_list[0]#Energy modification coefficient under DBE
gama_MCE  = float(gama_list[1].item())#Energy modification coefficient under MCE
# print('γDBE2=',gama_DBE,'γMCE4=',gama_MCE)
g =9.8 ##acceleration of gravity (N/kg)
if T<=T0:
    Sa_MCE = Sms*(0.4+0.6*(T/T0))
    Sa_DBE = (2/3)*Sa_MCE
if T0<T<=Ts:
    Sa_MCE = Sms
    Sa_DBE = (2/3)*Sa_MCE
if Ts<T<=4.0:
    Sa_MCE = Sm1/T
    Sa_DBE = (2/3)*Sa_MCE
if T>4:
    Sa_MCE = Sm1*4/(T**2)
    Sa_DBE = (2/3)*Sa_MCE
#Height of each subassemblage models (ST1 & ST2)
l=4600 ## mm
l3=1250## mm
h1 = 1400+1333+1400##The height of the first layer of subassemblage models
h2 = 4200##The height of the second layer of subassemblage models
h3 = 4200##The height of the third layer of subassemblage models
h4 = 4200##The height of the fourth layer of subassemblage models
h5 = 4200##The height of the fifth layer of subassemblage models
h6 = 4200-1400##The height of the sixth layer of subassemblage models
H_t = 26400 ###Total height of structure
N_damp = 20 ##Number of dampers per story
##Calculate the ratios of inelastic energy
eta_MCE = (2*(1-a1_)*(miu-1)+(a2_-a1_)*(ξ2_-ξ1_)*(2*miu-ξ1_-ξ2_)+(a3_-a1_)*(miu-ξ2_)**2)/(2*miu-1+a1_*(ξ1_-1)*(2*miu-ξ1_-1)+a2_*(ξ2_-ξ1_)*(2*miu-ξ1_-ξ2_)+a3_*(miu-ξ2_)**2)  ##塑性能量比
###Seismic weight and height to base for each story
W1 = 6348.0*10**3 #kN
W2 = 6348.0*10**3 #kN
W3 = 6348.0*10**3 #kN
W4 = 6348.0*10**3 #kN
W5 = 6348.0*10**3 #kN
W6 = 5290.0*10**3 #kN
H1 = 5400  #mm
H2 = 9600  #mm
H3 = 13800 #mm
H4 = 18000 #mm
H5 = 22200 #mm
H6 = 26400 #mm
##vertical distribution factor
beta1 = ((W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6)/(W6*H6))**(0.75*T**(-0.2))#F1
beta2 = ((W2*H2+W3*H3+W4*H4+W5*H5+W6*H6)/(W6*H6))**(0.75*T**(-0.2))#F2
beta3 = ((W3*H3+W4*H4+W5*H5+W6*H6)/(W6*H6))**(0.75*T**(-0.2))#F3
beta4 = ((W4*H4+W5*H5+W6*H6)/(W6*H6))**(0.75*T**(-0.2))#F4
beta5 = ((W5*H5+W6*H6)/(W6*H6))**(0.75*T**(-0.2))#F5
beta6 = ((W6*H6)/(W6*H6))**(0.75*T**(-0.2))#F6
Cv1 = (beta1-beta2)*(((W6*H6)/(W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6))**(0.75*T**(-0.2)))#F1
Cv2 = (beta2-beta3)*(((W6*H6)/(W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6))**(0.75*T**(-0.2)))#F2
Cv3 = (beta3-beta4)*(((W6*H6)/(W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6))**(0.75*T**(-0.2)))#F3
Cv4 = (beta4-beta5)*(((W6*H6)/(W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6))**(0.75*T**(-0.2)))#F4
Cv5 = (beta5-beta6)*(((W6*H6)/(W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6))**(0.75*T**(-0.2)))#F5
Cv6 = (beta6-0)*(((W6*H6)/(W1*H1+W2*H2+W3*H3+W4*H4+W5*H5+W6*H6))**(0.75*T**(-0.2)))    #F6
Ee= 1000*M*T**2*Sa_MCE**2*g**2/(8*(math.pi)**2)#Elastic energy  N·mm  N·mm
Ep = Ee*gama_MCE*eta_MCE ##the inelastic input energy
#Calculate base shear corresponding to DBE and MCE, take the maximum value as the design base shear
b_MCE = (8*math.pi**2*((theta4_/100-theta1_/100)+((a1_*(theta2_/100-theta1_/100)*(2*theta4_/100-theta1_/100-theta2_/100)+a2_*(theta3_/100-theta2_/100)*(2*theta4_/100-theta3_/100-theta2_/100)+a3_*(theta4_/100-theta3_/100)**2))/(2*theta1_/100))*(Cv1*H1+Cv2*H2+Cv3*H3+Cv4*H4+Cv5*H5+Cv6*H6))/(T**2)
Vy1_MCE=M/1000*(((-b_MCE+math.sqrt(b_MCE**2+4*gama_MCE*(Sa_MCE*g*1000)**2)))/2)
b_DBE = (8*math.pi**2*((theta2_-theta1_)/100+a1_*(theta2_-theta1_)**2/(2*theta1_*100))*(Cv1*H1+Cv2*H2+Cv3*H3+Cv4*H4+Cv5*H5+Cv6*H6))/(T**2)
Vy1_DBE=M/1000*(((-b_DBE+math.sqrt(b_DBE**2+4*gama_DBE*(Sa_DBE*g*1000)**2)))/2)
V1_max = max(Vy1_MCE,Vy1_DBE)
print('V1max=',V1_max)
##Design base shear
Vy1 = V1_max#v1
K = Vy1*100/(theta1_*H_t)
Vy2 = Vy1+((theta2_-theta1_)*H_t/100)*a1_*K#V2
Vy3 = Vy2+((theta3_-theta2_)*H_t/100)*a2_*K#V3
Vy4 = Vy3+((theta4_-theta3_)*H_t/100)*a3_*K#V4
####The first stage
###Lateral force distribution at each story
Fy1_1 = Cv1*Vy1  ##F1     N
Fy1_2 = Cv2*Vy1  ##F2     N
Fy1_3 = Cv3*Vy1  ##F3     N
Fy1_4 = Cv4*Vy1  ##F4     N
Fy1_5 = Cv5*Vy1  ##F5     N
Fy1_6 = Cv6*Vy1  ##F6     N
##Design story shear at each level
Vy1_1 = Fy1_1+Fy1_2+Fy1_3+Fy1_4+Fy1_5+Fy1_6  #F1
Vy1_2 = Fy1_2+Fy1_3+Fy1_4+Fy1_5+Fy1_6 #F2
Vy1_3 = Fy1_3+Fy1_4+Fy1_5+Fy1_6 #F3
Vy1_4 = Fy1_4+Fy1_5+Fy1_6 #F4
Vy1_5 = Fy1_5+Fy1_6 #F5
Vy1_6 = Fy1_6 #F6                       

#### Damper moment values at each story -Static equilibrium
My1_1 = Vy1_1*h1*l3/(l*N_damp)  ##N·mm  F1  
My1_2 = Vy1_2*h2*l3/(l*N_damp)  ##N·mm  F2  
My1_3 = Vy1_3*h3*l3/(l*N_damp)  ##N·mm  F3  
My1_4 = Vy1_4*h4*l3/(l*N_damp)  ##N·mm  F4  
My1_5 = Vy1_5*h5*l3/(l*N_damp)  ##N·mm  F5  
My1_6 = Vy1_6*h6*l3/(l*N_damp)  ##N·mm  F6  
####The Second stage
###Lateral force distribution at each story
Vy2_1 = Vy1_1  # N
Vy2_2 = Vy1_2  # N
Vy2_3 = Vy1_3  # N
Vy2_4 = Vy1_4  # N
Vy2_5 = Vy1_5  # N
Vy2_6 = Vy1_6  # N
#### Damper moment values at each story -Static equilibrium
My2_1 = My1_1  ##N·mm
My2_2 = My1_2  ##N·mm
My2_3 = My1_3  ##N·mm
My2_4 = My1_4  ##N·mm
My2_5 = My1_5  ##N·mm
My2_6 = My1_6  ##N·mm
####The third stage
###Lateral force distribution at each story
Fy3_1 = Cv1*Vy3  ## N
Fy3_2 = Cv2*Vy3  ## N
Fy3_3 = Cv3*Vy3  ## N
Fy3_4 = Cv4*Vy3  ## N
Fy3_5 = Cv5*Vy3  ## N
Fy3_6 = Cv6*Vy3  ## N
###Lateral force distribution at each story
Vy3_1 = Fy3_1+Fy3_2+Fy3_3+Fy3_4+Fy3_5+Fy3_6 #F1
Vy3_2 = Fy3_2+Fy3_3+Fy3_4+Fy3_5+Fy3_6  #F2
Vy3_3 = Fy3_3+Fy3_4+Fy3_5+Fy3_6 #F3
Vy3_4 = Fy3_4+Fy3_5+Fy3_6 #F4
Vy3_5 = Fy3_5+Fy3_6 #F5
Vy3_6 = Fy3_6 #F6
#### Damper moment values at each story -Static equilibrium
My3_1 = Vy3_1*h1*l3/(l*N_damp)  ##N·mm F1
My3_2 = Vy3_2*h2*l3/(l*N_damp)  ##N·mm F2
My3_3 = Vy3_3*h3*l3/(l*N_damp)  ##N·mm F3
My3_4 = Vy3_4*h4*l3/(l*N_damp)  ##N·mm F4
My3_5 = Vy3_5*h5*l3/(l*N_damp)  ##N·mm F5
My3_6 = Vy3_6*h6*l3/(l*N_damp)  ##N·mm F6
####The fourth stage
###Lateral force distribution at each story
Fy4_1 = Cv1*Vy4  ##F1     N
Fy4_2 = Cv2*Vy4  ##F2     N
Fy4_3 = Cv3*Vy4  ##F3     N
Fy4_4 = Cv4*Vy4  ##F4     N
Fy4_5 = Cv5*Vy4  ##F5     N
Fy4_6 = Cv6*Vy4  ##F6     N
###Lateral force distribution at each story
Vy4_1 = Fy4_1+Fy4_2+Fy4_3+Fy4_4+Fy4_5+Fy4_6  #F1
Vy4_2 = Fy4_2+Fy4_3+Fy4_4+Fy4_5+Fy4_6  #F2
Vy4_3 = Fy4_3+Fy4_4+Fy4_5+Fy4_6  #F3
Vy4_4 = Fy4_4+Fy4_5+Fy4_6 #F4
Vy4_5 = Fy4_5+Fy4_6 #F5
Vy4_6 = Fy4_6 #F6
#### Damper moment values at each story -Static equilibrium
My4_1 = Vy4_1*h1*l3/(l*N_damp)  ##N·mm F1
My4_2 = Vy4_2*h2*l3/(l*N_damp)  ##N·mm F2
My4_3 = Vy4_3*h3*l3/(l*N_damp)  ##N·mm F3
My4_4 = Vy4_4*h4*l3/(l*N_damp)  ##N·mm F4
My4_5 = Vy4_5*h5*l3/(l*N_damp)  ##N·mm F5
My4_6 = Vy4_6*h6*l3/(l*N_damp)  ##N·mm F6
#########################################
##Plastic energy distribution per story
#########################################
Ep1 = ((beta1)/(beta1+beta2+beta3+beta4+beta5+beta6))*Ep  ##N·mm F1
Ep2 = ((beta2)/(beta1+beta2+beta3+beta4+beta5+beta6))*Ep  ##N·mm F2
Ep3 = ((beta3)/(beta1+beta2+beta3+beta4+beta5+beta6))*Ep  ##N·mm F3
Ep4 = ((beta4)/(beta1+beta2+beta3+beta4+beta5+beta6))*Ep  ##N·mm F4
Ep5 = ((beta5)/(beta1+beta2+beta3+beta4+beta5+beta6))*Ep  ##N·mm F5
Ep6 = ((beta6)/(beta1+beta2+beta3+beta4+beta5+beta6))*Ep  ##N·mm F6
####Plastic energy provided by dampers at each story
##Calculate damper rotation angle   -Assume uniform inter-story drift ratio along height
ro1 = l*theta1_/(l3*100)  ##rad
ro2 = l*theta2_/(l3*100)
ro3 = l*theta3_/(l3*100)
ro4 = l*theta4_/(l3*100)
#The inelastic energy dissipation capacities
Ep1_damp = N_damp*(My1_1*(0.5*ro2+0.5*ro3-ro1)+My3_1*(0.5*ro4-0.5*ro2)+My4_1*(0.5*ro4-0.5*ro3))  ##N·mm·rad
Ep2_damp = N_damp*(My1_2*(0.5*ro2+0.5*ro3-ro1)+My3_2*(0.5*ro4-0.5*ro2)+My4_2*(0.5*ro4-0.5*ro3))
Ep3_damp = N_damp*(My1_3*(0.5*ro2+0.5*ro3-ro1)+My3_3*(0.5*ro4-0.5*ro2)+My4_3*(0.5*ro4-0.5*ro3))
Ep4_damp = N_damp*(My1_4*(0.5*ro2+0.5*ro3-ro1)+My3_4*(0.5*ro4-0.5*ro2)+My4_4*(0.5*ro4-0.5*ro3))
Ep5_damp = N_damp*(My1_5*(0.5*ro2+0.5*ro3-ro1)+My3_5*(0.5*ro4-0.5*ro2)+My4_5*(0.5*ro4-0.5*ro3))
Ep6_damp = N_damp*(My1_6*(0.5*ro2+0.5*ro3-ro1)+My3_6*(0.5*ro4-0.5*ro2)+My4_6*(0.5*ro4-0.5*ro3))
#Initial damper plastic energy dissipation capacity is insufficient. Increase damper capacity to meet plastic energy demand.
f1_damp = Ep1/Ep1_damp
f2_damp = Ep2/Ep2_damp
f3_damp = Ep3/Ep3_damp
f4_damp = Ep4/Ep4_damp
f5_damp = Ep5/Ep5_damp
f6_damp = Ep6/Ep6_damp
if f1_damp<1:
    f1_damp =1
if f2_damp<1:
    f2_damp =1
if f3_damp<1:
    f3_damp =1
if f4_damp<1:
    f4_damp =1
if f5_damp<1:
    f5_damp =1
if f6_damp<1:
    f6_damp =1
Ep1_damp = f1_damp*N_damp*(My1_1*(0.5*ro2+0.5*ro3-ro1)+My3_1*(0.5*ro4-0.5*ro2)+My4_1*(0.5*ro4-0.5*ro3))  ##N·mm·rad
Ep2_damp = f2_damp*N_damp*(My1_2*(0.5*ro2+0.5*ro3-ro1)+My3_2*(0.5*ro4-0.5*ro2)+My4_2*(0.5*ro4-0.5*ro3))
Ep3_damp = f3_damp*N_damp*(My1_3*(0.5*ro2+0.5*ro3-ro1)+My3_3*(0.5*ro4-0.5*ro2)+My4_3*(0.5*ro4-0.5*ro3))
Ep4_damp = f4_damp*N_damp*(My1_4*(0.5*ro2+0.5*ro3-ro1)+My3_4*(0.5*ro4-0.5*ro2)+My4_4*(0.5*ro4-0.5*ro3))
Ep5_damp = f5_damp*N_damp*(My1_5*(0.5*ro2+0.5*ro3-ro1)+My3_5*(0.5*ro4-0.5*ro2)+My4_5*(0.5*ro4-0.5*ro3))
Ep6_damp = f6_damp*N_damp*(My1_6*(0.5*ro2+0.5*ro3-ro1)+My3_6*(0.5*ro4-0.5*ro2)+My4_6*(0.5*ro4-0.5*ro3))
# print(f1_damp,f2_damp,f3_damp,f4_damp,f5_damp,f6_damp)

###Damper stiffness, capacity, slip angle
s = ro2-ro1  ##Damper friction segment slip angle
###F1
##the design moment Mji of the damper
My1_1_d = My1_1*f1_damp  #first stage N·mm
My2_1_d = My1_1_d        #second stage
My3_1_d = My3_1*f1_damp  #third stage
My4_1_d = My4_1*f1_damp  #fourth stage
###the design stiffness of the damper
k1_1 = My1_1_d/ro1  ##N·mm/rad
k2_1 = 0
k3_1 = (My3_1_d-My1_1_d)/(ro3-ro2)
k4_1 = (My4_1_d-My3_1_d)/(ro4-ro3)
###F2
##the design moment Mji of the damper
My1_2_d = My1_2*f2_damp #first stage N·mm 
My2_2_d = My1_2_d       #second stage
My3_2_d = My3_2*f2_damp #third stage
My4_2_d = My4_2*f2_damp #fourth stage
###the design stiffness of the damper
k1_2 = My1_2_d/ro1  ##N·mm/rad
k2_2 = 0
k3_2 = (My3_2_d-My1_2_d)/(ro3-ro2)
k4_2 = (My4_2_d-My3_2_d)/(ro4-ro3)
###F3
##the design moment Mji of the damper
My1_3_d = My1_3*f3_damp  #first stage N·mm 
My2_3_d = My1_3_d        #second stage
My3_3_d = My3_3*f3_damp  #third stage
My4_3_d = My4_3*f3_damp  #fourth stage
###the design stiffness of the damper
k1_3 = My1_3_d/ro1  ##N·mm/rad
k2_3 = 0
k3_3 = (My3_3_d-My1_3_d)/(ro3-ro2)
k4_3 = (My4_3_d-My3_3_d)/(ro4-ro3)
###F4
##the design moment Mji of the damper
My1_4_d = My1_4*f4_damp  #first stage N·mm 
My2_4_d = My1_4_d        #second stage
My3_4_d = My3_4*f4_damp  #third stage
My4_4_d = My4_4*f4_damp  #fourth stage
###the design stiffness of the damper
k1_4 = My1_4_d/ro1  ##N·mm/rad
k2_4 = 0
k3_4 = (My3_4_d-My1_4_d)/(ro3-ro2)
k4_4 = (My4_4_d-My3_4_d)/(ro4-ro3)
###F5
My1_5_d = My1_5*f5_damp ##N·mm
My2_5_d = My1_5_d  
My3_5_d = My3_5*f5_damp 
My4_5_d = My4_5*f5_damp 
k1_5 = My1_5_d/ro1  ##N·mm/rad
k2_5 = 0
k3_5 = (My3_5_d-My1_5_d)/(ro3-ro2)
k4_5 = (My4_5_d-My3_5_d)/(ro4-ro3)
###F6
My1_6_d = My1_6*f6_damp   ##N·mm
My2_6_d = My1_6_d  
My3_6_d = My3_6*f6_damp 
My4_6_d = My4_6*f6_damp 
k1_6 = My1_6_d/ro1  ##N·mm/rad
k2_6 = 0
k3_6 = (My3_6_d-My1_6_d)/(ro3-ro2)
k4_6 = (My4_6_d-My3_6_d)/(ro4-ro3)
# print('F1_damp=',{'Kf':k1_1,'Kr':k3_1 ,'Kp':k4_1,'Mf':My1_1_d,'Mr':My3_1_d,'s':s})
# print('F2_damp=',{'Kf':k1_2,'Kr':k3_2 ,'Kp':k4_2,'Mf':My1_2_d,'Mr':My3_2_d,'s':s})
# print('F3_damp=',{'Kf':k1_3,'Kr':k3_3 ,'Kp':k4_3,'Mf':My1_3_d,'Mr':My3_3_d,'s':s})
# print('F4_damp=',{'Kf':k1_4,'Kr':k3_4 ,'Kp':k4_4,'Mf':My1_4_d,'Mr':My3_4_d,'s':s})
# print('F5_damp=',{'Kf':k1_5,'Kr':k3_5 ,'Kp':k4_5,'Mf':My1_5_d,'Mr':My3_5_d,'s':s})
# print('F6_damp=',{'Kf':k1_6,'Kr':k3_6 ,'Kp':k4_6,'Mf':My1_6_d,'Mr':My3_6_d,'s':s})
# print('Mp_list =',[My4_1_d,My4_2_d,My4_3_d,My4_4_d,My4_5_d,My4_6_d])
##The required plastic moment capacity of the column base
Vy = 4*3.14159**2*M*theta1_*(Cv1*H1+Cv2*H2+Cv3*H3+Cv4*H4+Cv5*H5+Cv6*H6)*10**(-3)/(T**2)
Mpc =(Vy4/5)*1.1*H1/4
Wp1=Mpc/420
Wp2=2*Mpc/420
# print('Mpc=',Mpc)
print('Wp1=',Wp1)
print('Wp2=',Wp2)
##Select the section of the side column and the middle column.
Mpc1 = 5229238.824*420 #the side column
Mpc2 = 11002407.679000003*420 #the middle column
Mpc_all = 2*Mpc1+4*Mpc2 #the sum of the plastic moment capacity of the column bases
Mdampc = My4_1_d+My4_2_d+My4_3_d+My4_4_d+My4_5_d+My4_6_d #the sum of the the TSRD moment at the 4th segment
Vp4 = (Mpc_all+20*(l/l3)*Mdampc)/(Cv1*H1+Cv2*H2+Cv3*H3+Cv4*H4+Cv5*H5+Cv6*H6)
print('V4=',Vp4/1000)

############################
##According to the theoretical formula of subassemblage, check the actual structural parameters.
############################
###The height of each section of substructure [Floor1,Floor2,Floor3,Floor4,Floor5,Floor6]
h_list  = [2800+1333,4200,4200,4200,4200,4200-1400] #total height
h1_list = [1400,1400,1400,1400,1400,0]   #h1
h2_list = [1400,1400,1400,1400,1400,1400]#h2
h3_list = [1333,1400,1400,1400,1400,1400]#h3
###Section modulus and area of truss members and columns [Floor1,Floor2,Floor3,Floor4,Floor5,Floor6]
I1_list_ST1 = [401717574.0, 401717574.0, 368590885.3333333, 368590885.3333333, 347257552.0, 347257552.0]
A1_list_ST1 = [18402.0, 18402.0, 16444.0, 16444.0, 14844.0, 14844.0]
I3_list_ST1 = [401717574.0, 401717574.0, 368590885.3333333, 368590885.3333333, 347257552.0, 347257552.0]
A3_list_ST1 = [18402.0, 18402.0, 16444.0, 16444.0, 14844.0, 14844.0]
A5_list_ST1 = [18402.0, 18402.0, 16444.0, 16444.0, 14844.0, 14844.0]
I7_list_ST1 = I8_list_ST1 = [889307450.5031998, 889307450.5031998, 783046051.9258667, 783046051.9258667, 706655767.4157333, 706655767.4157333]
A7_list_ST1 = A8_list_ST1 = [33267.24, 33267.24, 29898.44, 29898.44, 27356.68, 27356.68]
I7_list_ST2 = I8_list_ST2 = [2040518269.1624668, 2040518269.1624668, 1796871226.9091997, 1796871226.9091997, 1591076525.8274672, 1591076525.8274672]
A7_list_ST2 = A8_list_ST2 = [64701.86, 64701.86, 58760.759999999995, 58760.759999999995, 53516.240000000005, 53516.240000000005]
I2_list_ST1 =I1_list_ST1
A2_list_ST1 =A1_list_ST1
A4_list_ST1 =A3_list_ST1
theta1 = 0.28
###Initial input of damper strength and stiffness information
F1_damp= {'Kf': k1_1, 'Kr': k3_1 , 'Kp': k4_1, 'Mf': My1_1_d, 'Mr': My3_1_d, 's': s}
F2_damp= {'Kf': k1_2, 'Kr': k3_2 , 'Kp': k4_2, 'Mf': My1_2_d, 'Mr': My3_2_d, 's': s}
F3_damp= {'Kf': k1_3, 'Kr': k3_3 , 'Kp': k4_3, 'Mf': My1_3_d, 'Mr': My3_3_d, 's': s}
F4_damp= {'Kf': k1_4, 'Kr': k3_4 , 'Kp': k4_4, 'Mf': My1_4_d, 'Mr': My3_4_d, 's': s}
F5_damp= {'Kf': k1_5, 'Kr': k3_5 , 'Kp': k4_5, 'Mf': My1_5_d, 'Mr': My3_5_d, 's': s}
F6_damp= {'Kf': k1_6, 'Kr': k3_6 , 'Kp': k4_6, 'Mf': My1_6_d, 'Mr': My3_6_d, 's': s}
Mp_list = [My4_1_d,My4_2_d,My4_3_d,My4_4_d,My4_5_d,My4_6_d]
theta_list = []
K_list =[]
beishu_list  =[1,1,1,1,1,1]
##Adjust damper stiffness to meet displacement target requirements
beishu_list1=[2.0328411045644943, 2.0211101823954127, 2.0490442469868584, 1.9171704056913625, 1.7981730938900191, 1.4712680434502556]
beishu_list2 =[1.488472207904802, 1.4798004440522883, 1.500507542341603, 1.4046011658318063, 1.3225873671070274, 1.1280742598282456]
beishu_list=np.multiply(beishu_list1, beishu_list2)
beishu_list3 =[1.3155069340725565, 1.307918314581937, 1.3260923423048745, 1.2436194276777575, 1.177257554314233, 1.0453951891758353]
beishu_list=np.multiply(beishu_list3, beishu_list)
beishu_kp =0.363
for Snum in range(1,7): 
    ######ST1 formula
    h =h_list[Snum-1]
    h1=h1_list[Snum-1]
    h2=h2_list[Snum-1]
    h3=h3_list[Snum-1]
    ##Chord member lengths of truss beam segments
    l =4600
    l1=1675
    l2=1675
    l3=1250
    I1=I1_list_ST1[Snum-1]
    I2=I2_list_ST1[Snum-1]
    I3=I3_list_ST1[Snum-1]
    I7=I7_list_ST1[Snum-1]
    I8=I8_list_ST1[Snum-1]
    A1=A1_list_ST1[Snum-1]
    A3=A3_list_ST1[Snum-1]
    A4=A4_list_ST1[Snum-1]
    A5=A5_list_ST1[Snum-1]
    A7=A7_list_ST1[Snum-1]
    A8=A8_list_ST1[Snum-1]
    E=2.0*10**5    #Elastic modulus
    u = 0.3        #Poisson's ratio
    G = E/(2*(1+u))#Shear modulus
    #the strength and stiffness of dampers
    F_damp_list =[F1_damp,F2_damp,F3_damp,F4_damp,F5_damp,F6_damp] 
    F_damp = F_damp_list[Snum-1]
    k_damp_e = F_damp['Kf']*beishu_list[Snum-1]
    k_damp_y = F_damp['Kr']*beishu_list[Snum-1]
    k_damp_p = F_damp['Kp']*beishu_list[Snum-1]*beishu_kp
    Me = F_damp['Mf']
    My = F_damp['Mr']
    Mp =Mp_list[Snum-1]
    s = F_damp['s']
    ###calculation coefficient
    C1 = h2/(3*E*I2)+l2/(3*E*I3)+l2/(E*A3*h2**2)+(l2**2+h2**2)**(3/2)/(E*A5*l2**2*h2**2)
    C2=h2/(6*E*I2)-l2/(E*A3*h2**2)-(l2**2+h2**2)**(3/2)/(E*A5*l2**2*h2**2)
    C3=h*(l2+l3)*(h2**2+l2**2)**(3/2)/(E*A5*l*l2**2*h2**2)
    C4=l2*l3/(3*E*I3)
    X1 = h/(2*l)
    X2 = (C3+X1*C4)/(C2-C1)
    X3=-X2
    F1 = math.sqrt(l2**2+h2**2)/(l2*h2)
    F2 = (h*(l2+l3)*math.sqrt(l2**2+h2**2))/(l*l2*h2)
    F3 = h*(l2+l3)/(l*h2)
    F4 = h*math.sqrt(h2**2+l1**2)/(l*h2)
    ###Displacement - horizontal displacement at column top under unit load
    deta_M = (h**2*l3**3)/(6*E*I1*l**2)+(l2*l3*h)*(h*l3/(2*l)-X3)/(3*E*I3*l)+(h1**2*h2/3-h1*h2*h3/3+h2*h3**2/3+h3**3/3)/(E*I7)+h1**3/(3*E*I8)
    deta_Q = l1*h**2/(4*G*A1*l**2)+(h3+(h1+h3)**2/h2)/(G*A7)+h1/(G*A8)
    deta_N = F3**2*l1/(E*A3)+l1*(F3**2+((h1+h3)/h2+h/l)**2)/(E*A4)+((F2+F1*X2-F1*X3)*F2*math.sqrt(h2**2+l2**2)+F4**2*math.sqrt(h2**2+l1**2))/(E*A5)
    deta_damp_e = h**2*l3**2/(2*l**2*k_damp_e)
    deta_damp_y = h**2*l3**2/(2*l**2*k_damp_y)
    deta_damp_p = h**2*l3**2/(2*l**2*k_damp_p)
    ###subassemblage stiffness
    K_STMF1_e = h/(deta_M+deta_Q+deta_N+deta_damp_e)
    K_STMF1_y = h/(deta_M+deta_Q+deta_N+deta_damp_y)
    K_STMF1_p = h/(deta_M+deta_Q+deta_N+deta_damp_p)
    ####subassemblage slip drift ratio
    theta1_d = s*l3/l
    ###subassemblage column shear
    V_STMF1_e=2*l*Me/(h*l3)
    V_STMF1_y=2*l*My/(h*l3)
    V_STMF1_p=2*l*Mp/(h*l3)
    theta1_e = V_STMF1_e/K_STMF1_e
    theta1_s = theta1_e+theta1_d
    theta1_y = theta1_s+(V_STMF1_y-V_STMF1_e)/K_STMF1_y
    theta1_p = theta1_y+(V_STMF1_p-V_STMF1_y)/K_STMF1_p
    ######ST2 formula
    I7=I7_list_ST2[Snum-1]
    I8=I8_list_ST2[Snum-1]
    A7=A7_list_ST2[Snum-1]
    A8=A8_list_ST2[Snum-1]
    ###calculation coefficient
    C1 = h2/(3*E*I2)+l2/(3*E*I3)+l2/(E*A3*h2**2)+(l2**2+h2**2)**(3/2)/(E*A5*l2**2*h2**2)
    C2 =h2/(6*E*I2)-l2/(E*A3*h2**2)-(l2**2+h2**2)**(3/2)/(E*A5*l2**2*h2**2)
    C3 =h*(l2+l3)*(h2**2+l2**2)**(3/2)/(E*A5*l*l2**2*h2**2)
    C4 =l2*l3/(3*E*I3)
    X1 = h/(2*l)
    X2 = (C3+X1*C4)/(C2-C1)
    X3 =-X2
    F1 = math.sqrt(l2**2+h2**2)/(l2*h2)
    F2 = (h*(l2+l3)*math.sqrt(l2**2+h2**2))/(l*l2*h2)
    F3 = h*(l2+l3)/(l*h2)
    F4 = h*math.sqrt(h2**2+l1**2)/(l*h2)
    ###Displacement - horizontal displacement at column top under unit load
    deta_M=(h**2*l3**3)/(12*E*I1*l**2)+(l2*l3*h)*(h*l3/(2*l)-X3)/(6*E*I3*l)+(h1**2*h2/3-h1*h2*h3/3+h2*h3**2/3+h3**3/3)/(E*I7)+h1**3/(3*E*I8)
    deta_Q=l1*h**2/(8*G*A1*l**2)+(h3+(h1+h3)**2/h2)/(G*A7)+h1/(G*A8)
    deta_N = F3**2*l1/(2*E*A3)+l1*(F3**2+((h1+h3)/h2+h/l)**2)/(2*E*A4)+((F2+F1*X2-F1*X3)*F2*math.sqrt(h2**2+l2**2)+F4**2*math.sqrt(h2**2+l1**2))/(2*E*A5)
    deta_damp_e = h**2*l3**2/(4*l**2*k_damp_e)
    deta_damp_y = h**2*l3**2/(4*l**2*k_damp_y)
    deta_damp_p = h**2*l3**2/(4*l**2*k_damp_p)
    ###subassemblage stiffness
    K_STMF2_e = h/(deta_M+deta_Q+deta_N+deta_damp_e)
    K_STMF2_y = h/(deta_M+deta_Q+deta_N+deta_damp_y)
    K_STMF2_p = h/(deta_M+deta_Q+deta_N+deta_damp_p)
    ####subassemblage slip drift ratio
    theta2_d = s*l3/l
    ###subassemblage column shear
    V_STMF2_e=4*l*Me/(h*l3)
    V_STMF2_y=4*l*My/(h*l3)
    V_STMF2_p=4*l*Mp/(h*l3)
    theta2_e = V_STMF2_e/K_STMF2_e
    theta2_s = theta2_e+theta2_d
    theta2_y = theta2_s+(V_STMF2_y-V_STMF2_e)/K_STMF2_y
    theta2_p = theta2_y+(V_STMF2_p-V_STMF2_y)/K_STMF2_p
    ##The number of subassemblages in each layer
    N_ST1 = 2
    N_ST2 = 4
    ##Middle layer stiffness
    K_mid_e = N_ST1*K_STMF1_e+N_ST2*K_STMF2_e
    K_mid_y = N_ST1*K_STMF1_y+N_ST2*K_STMF2_y
    K_mid_p = N_ST1*K_STMF1_p+N_ST2*K_STMF2_p
    theta_mid_e = max(theta1_e,theta2_e)
    theta_mid_s = max(theta1_s,theta2_s)
    theta_mid_y = max(theta1_y,theta2_y)
    theta_mid_p = max(theta1_p,theta2_p)

    V_mid_e = K_mid_e*theta_mid_e
    V_mid_y = K_mid_y*(theta_mid_y-theta_mid_s)+V_mid_e
    V_mid_p = K_mid_p*(theta_mid_p-theta_mid_y)+V_mid_y
    # print('x_theory = ',[0,theta_mid_e*100,theta_mid_s*100,theta_mid_y*100,theta_mid_p*100])
    # print('y_theory =' ,[0,V_mid_e/1000,V_mid_e/1000,V_mid_y/1000,V_mid_p/1000])
    # print(K_mid_e,K_mid_y,K_mid_p)
    K_list.append(theta_mid_e*100/theta1)
    theta_list.append([0,theta_mid_e*100,theta_mid_s*100,theta_mid_y*100,theta_mid_p*100])

# print(K_list)
theta_0=(theta_list[0][0]+theta_list[1][0]+theta_list[2][0]+theta_list[3][0]+theta_list[4][0]+theta_list[5][0])/6
theta_1=(theta_list[0][1]+theta_list[1][1]+theta_list[2][1]+theta_list[3][1]+theta_list[4][1]+theta_list[5][1])/6
theta_2=(theta_list[0][2]+theta_list[1][2]+theta_list[2][2]+theta_list[3][2]+theta_list[4][2]+theta_list[5][2])/6
theta_3=(theta_list[0][3]+theta_list[1][3]+theta_list[2][3]+theta_list[3][3]+theta_list[4][3]+theta_list[5][3])/6
theta_4=(theta_list[0][4]+theta_list[1][4]+theta_list[2][4]+theta_list[3][4]+theta_list[4][4]+theta_list[5][4])/6
###Output theoretical drift ratio of structure [0, θ1, θ2, θ3, θ4]
print([theta_0,theta_1,theta_2,theta_3,theta_4])
###Output final damper design results  N mm
######
print('Dampers design information')
print('F1_damp=',{'M1':My1_1_d,'M3':My3_1_d,'M4':My4_1_d,'K1':k1_1*beishu_list[0],'K3':k3_1*beishu_list[0] ,'K4':k4_1*beishu_list[0]*beishu_kp,'S':s})
print('F2_damp=',{'M1':My1_2_d,'M3':My3_2_d,'M4':My4_2_d,'K1':k1_2*beishu_list[1],'K3':k3_2*beishu_list[1] ,'K4':k4_2*beishu_list[1]*beishu_kp,'S':s})
print('F3_damp=',{'M1':My1_3_d,'M3':My3_3_d,'M4':My4_3_d,'K1':k1_3*beishu_list[2],'K3':k3_3*beishu_list[2] ,'K4':k4_3*beishu_list[2]*beishu_kp,'S':s})
print('F4_damp=',{'M1':My1_4_d,'M3':My3_4_d,'M4':My4_4_d,'K1':k1_4*beishu_list[3],'K3':k3_4*beishu_list[3] ,'K4':k4_4*beishu_list[3]*beishu_kp,'S':s})
print('F5_damp=',{'M1':My1_5_d,'M3':My3_5_d,'M4':My4_5_d,'K1':k1_5*beishu_list[4],'K3':k3_5*beishu_list[4] ,'K4':k4_5*beishu_list[4]*beishu_kp,'S':s})
print('F6_damp=',{'M1':My1_6_d,'M3':My3_6_d,'M4':My4_6_d,'K1':k1_6*beishu_list[5],'K3':k3_6*beishu_list[5] ,'K4':k4_6*beishu_list[5]*beishu_kp,'S':s})