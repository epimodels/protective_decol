#########################################################
# Dynamic Transmission Model of MRSA Cohorts.      	#
# Queue-based Steady State Populations             	#
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)   #
#########################################################

# Descriptive Information for PML File
Modelname: MRSA Cohort
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) Cohort 1 #
R1:
	N_u1 > N_c1
	rho_N * sigma * (P_c1+P_i1) * (N_u1 / (P_c1 + P_u1 + P_p1 +P_i1))

R2:
	N_c1 > N_u1
	N_c1 * iota_N

R3:
	N_c1 > N_u1
	N_c1 * tau_N * ((P_c1+P_i1) / (P_c1 + P_u1 + P_p1 +P_i1))
	

# Reactions Governing Movement of Nurses (N) Cohort 2 #
R4:
	N_u2 > N_c2
	rho_N * sigma * (P_c2+P_i2) * (N_u2 / (P_c2 + P_u2 +P_p2 + P_i2))

R5:
	N_c2 > N_u2
	N_c2 * iota_N

R6:
	N_c2 > N_u2
	N_c2 * tau_N * ((P_c2+P_i2) / (P_c2 + P_u2 + P_p2+P_i2))

# Reactions Governing Movement of Nurses (N) Cohort 3 #
R7:
	N_u3 > N_c3
	rho_N * sigma * (P_c3+P_i3) * (N_u3 / (P_c3 + P_u3 + P_p3+P_i3))

R8:
	N_c3 > N_u3
	N_c3 * iota_N

R9:
	N_c3 > N_u3
	N_c3 * tau_N * ((P_c3+P_i3) / (P_c3 + P_u3 + P_p3 +P_i3))

# Reactions Governing Movement of Nurses (N) Cohort 4 #
R10:
	N_u4 > N_c4
	rho_N * sigma * (P_c4+P_i4) * (N_u4 / (P_c4 + P_u4 + P_p4 +P_i4))

R11:
	N_c4 > N_u4
	N_c4 * iota_N

R12:
	N_c4 > N_u4
	N_c4 * tau_N * ((P_c4+P_i4) / (P_c4 + P_u4 + P_p4+P_i4))

# Reactions Governing Movement of Nurses (N) Cohort 5 #
R13:
	N_u5 > N_c5
	rho_N * sigma * (P_c5+P_i5) * (N_u5 / (P_c5 + P_u5 + P_p5+P_i5))

R14:
	N_c5 > N_u5
	N_c5 * iota_N

R15:
	N_c5 > N_u5
	N_c5 * tau_N * ((P_c5+P_i5) / (P_c5 + P_u5 + P_p5+P_i5))

# Reactions Governing Movement of Nurses (N) Cohort 6 #
R16:
	N_u6 > N_c6
	rho_N * sigma * (P_c6+P_i6) * (N_u6 / (P_c6 + P_u6 + P_p6+P_i6))

R17:
	N_c6 > N_u6
	N_c6 * iota_N

R18:
	N_c6 > N_u6
	N_c6 * tau_N * ((P_c6+P_i6) / (P_c6 + P_u6 + P_p6+P_i6))

########################

# Reactions Governing Movement of the Doctor 
R19:
	D_u > D_c
	rho_D * sigma * (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6+P_i1+P_i2+P_i3+P_i4+P_i5+P_i6) * (D_u / (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6 + P_u1+P_u2+P_u3+P_u4+P_u5+P_u6+ P_p1+P_p2+P_p3+P_p4+P_p5+P_p6+P_i1+P_i2+P_i3+P_i4+P_i5+P_i6))

R20:
	D_c > D_u
	D_c * iota_D

R21:
	D_c > D_u
	D_c * tau_D * ((P_c1+P_c2+P_c3+P_c4+P_c5+P_c6+P_i1+P_i2+P_i3+P_i4+P_i5+P_i6) / (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6 + P_u1+P_u2+P_u3+P_u4+P_u5+P_u6 +P_p1+P_p2+P_p3+P_p4+P_p5+P_p6 + P_p1+P_p2+P_p3+P_p4+P_p5+P_p6+P_i1+P_i2+P_i3+P_i4+P_i5+P_i6))
	

########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #
R22:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1/(N_u1 + N_c1))

R23:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c / (D_u + D_c))
	
R24:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R25:	
	P_u1 > P_c1
	theta * P_u1 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #
R26:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2/(N_u2 + N_c2))

R27:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c / (D_u + D_c))
	
R28:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R29:	
	P_u2 > P_c2
	theta * P_u2 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #
R30:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3/(N_u3 + N_c3))

R31:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c / (D_u + D_c))
	
R32:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R33:	
	P_u3 > P_c3
	theta * P_u3 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #
R34:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c4/(N_u4 + N_c4))

R35:
	P_u4 > P_c4 + Acquisition
	rho_D * psi * P_u4 * (D_c / (D_u + D_c))
	
R36:
	P_u4 > P_u4
	theta * P_u4 * (1-nu)

R37:	
	P_u4 > P_c4
	theta * P_u4 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #
R38:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c5/(N_u5 + N_c5))

R39:
	P_u5 > P_c5 + Acquisition
	rho_D * psi * P_u5 * (D_c / (D_u + D_c))
	
R40:
	P_u5 > P_u5
	theta * P_u5 * (1-nu)

R41:	
	P_u5 > P_c5
	theta * P_u5 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 6 #
R42:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c6/(N_u6 + N_c6))

R43:
	P_u6 > P_c6 + Acquisition
	rho_D * psi * P_u6 * (D_c / (D_u + D_c))
	
R44:
	P_u6 > P_u6
	theta * P_u6 * (1-nu)

R45:	
	P_u6 > P_c6
	theta * P_u6 * nu

########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #	
R46:
    P_c1 > P_u1
    mu*P_c1
    
R47:
	P_c1 > P_c1
	theta * P_c1 * nu * (1-deltaP)

R48:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #
R49:
    P_c2 > P_u2
    mu*P_c2
    

R50:
	P_c2 > P_c2
	theta * P_c2 * nu * (1-deltaP)

R51:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #
R52:
    P_c3 > P_u3
    mu * P_c3
    
R53:
	P_c3 > P_c3
	theta * P_c3 * nu * (1-deltaP)

R54:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #
R55:
    P_c4 > P_u4
    mu * P_c4
    
R56:
	P_c4 > P_c4
	theta * P_c4 * nu * (1-deltaP)

R57:
	P_c4 > P_u4
	theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #
R58:
    P_c5 > P_u5
    mu * P_c5
    
R59:
	P_c5 > P_c5
	theta * P_c5 * nu * (1-deltaP)

R60:
	P_c5 > P_u5
	theta * P_c5 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 6 #
R61:
    P_c6 > P_u6
    mu * P_c6
    
R62:
	P_c6 > P_c6
	theta * P_c6 * nu * (1-deltaP)

R63:
	P_c6 > P_u6
	theta * P_c6 * (1-nu)
	
# CHG-related Reactions
	
R64:
    P_c1 > P_p1
    epsilonE * epsilonB * P_c1
    
R65:
    P_u1 > P_p1
    epsilonE * epsilonB * P_u1
    
R66:
    P_p1 > P_u1
    epsilonU * P_p1
    
R67:
    P_c2 > P_p2
    epsilonE * epsilonB * P_c2
    
R68:
    P_u2 > P_p2
    epsilonE * epsilonB * P_u2
    
R69:
    P_p2 > P_u2
    epsilonU * P_p2
    
R70:
    P_c3 > P_p3
    epsilonE * epsilonB * P_c3
    
R71:
    P_u3 > P_p3
    epsilonE * epsilonB * P_u3
    
R72:
    P_p3 > P_u3
    epsilonU * P_p3
    
R73:
    P_c4 > P_p4
    epsilonE * epsilonB * P_c4
    
R74:
    P_u4 > P_p4
    epsilonE * epsilonB * P_u4
    
R75:
    P_p4 > P_u4
    epsilonU * P_p4
    
R76:
    P_c5 > P_p5
    epsilonE * epsilonB * P_c5
    
R77:
    P_u5 > P_p5
    epsilonE * epsilonB * P_u5
    
R78:
    P_p5 > P_u5
    epsilonU * P_p5
    
R79:
    P_c6 > P_p6
    epsilonE * epsilonB * P_c6
    
R80:
    P_u6 > P_p6
    epsilonE * epsilonB * P_u6
    
R81:
    P_p6 > P_u6
    epsilonU * P_p6
    
## Admission reactions for CHG
R82:
    P_p1 > P_u1
    mu * P_p1

R83:
	P_p1 > P_u1
	theta * P_p1 * (1-nu)
	
R84:
    P_p2 > P_u2
    mu * P_p2

R85:
	P_p2 > P_u2
	theta * P_p2 * (1-nu)
	
R86:
    P_p3 > P_u3
    mu * P_p3

R87:
	P_p3 > P_u3
	theta * P_p3 * (1-nu)
	
R88:
    P_p4 > P_u4
    mu * P_p4

R89:
	P_p4 > P_u4
	theta * P_p4 * (1-nu)
	
R90:
    P_p5 > P_u5
    mu * P_p5

R91:
	P_p5 > P_u5
	theta * P_p5 * (1-nu)
	
R92:
    P_p6 > P_u6
    mu * P_p6

R93:
	P_p6 > P_u6
	theta * P_p6 * (1-nu)
	
## Reactions involving infected individuals ##

## Reactions involving becoming infected ##

R94:
    P_c1 > P_i1 + Acquisition2
    deltaI * P_c1

R95:
    P_c2 > P_i2 + Acquisition2
    deltaI * P_c2

R96:
    P_c3 > P_i3 + Acquisition2
    deltaI * P_c3

R97:
    P_c4 > P_i4 + Acquisition2
    deltaI * P_c4

R98:
    P_c5 > P_i5 + Acquisition2
    deltaI * P_c5

R99:
    P_c6 > P_i6 + Acquisition2
    deltaI * P_c6

## Rate infected die ##

R100:
    P_i1 > P_u1 + Deaths
    deltaD * (1- nu)* P_i1
    
R101:
    P_i1 > P_c1 + Deaths
    deltaD * nu * P_i1
    
R102:
    P_i2 > P_u2 + Deaths
    deltaD * (1- nu)* P_i2

R103:
    P_i2 > P_c2 + Deaths
    deltaD * nu * P_i2

R104:
    P_i3 > P_u3 + Deaths
    deltaD * (1- nu)* P_i3

R105:
    P_i3 > P_c3 + Deaths
    deltaD * nu * P_i3

R106:
    P_i4 > P_u4 + Deaths
    deltaD * (1- nu)* P_i4

R107:
    P_i4 > P_c4 + Deaths
    deltaD * nu * P_i4

R108:
    P_i5 > P_u5 + Deaths
    deltaD * (1- nu)* P_i5

R109:
    P_i5 > P_c5 + Deaths
    deltaD * nu * P_i5

R110:
    P_i6 > P_u6 + Deaths
    deltaD * (1- nu)* P_i6

R111:
    P_i6 > P_c6 + Deaths
    deltaD * nu * P_i6
    
R112:
    P_i1 > P_u1
    deltaC * P_i1
    
R113:
    P_i2 > P_u2
    deltaC * P_i2
    
R114:
    P_i3 > P_u3
    deltaC * P_i3
    
R115:
    P_i4 > P_u4
    deltaC * P_i4
    
R116:
    P_i5 > P_u5
    deltaC * P_i5
    
R117:
    P_i6 > P_u6
    deltaC * P_i6
    
# admission with infection
R118:
	P_u1 > P_i1
	theta * P_u1 * nu * deltaP
    
R119:
	P_c1 > P_i1
	theta * P_c1 * nu * deltaP
	
R120:
	P_p1 > P_i1
	theta * P_p1 * nu * deltaP
	
R121:
	P_i1 > P_i1
	deltaD * P_i1 * nu * deltaP
	
R122:
	P_u2 > P_i2
	theta * P_u2 * nu * deltaP
    
R123:
	P_c2 > P_i2
	theta * P_c2 * nu * deltaP
	
R124:
	P_p2 > P_i2
	theta * P_p2 * nu * deltaP
	
R125:
	P_i2 > P_i2
	deltaD * P_i2 * nu * deltaP
	
R126:
	P_u3 > P_i3
	theta * P_u3 * nu * deltaP
    
R127:
	P_c3 > P_i3
	theta * P_c3 * nu * deltaP
	
R128:
	P_p3 > P_i3
	theta * P_p3 * nu * deltaP
	
R129:
	P_i3 > P_i3
	deltaD * P_i3 * nu * deltaP

R130:
	P_u4 > P_i4
	theta * P_u4 * nu * deltaP
    
R131:
	P_c4 > P_i4
	theta * P_c4 * nu * deltaP
	
R132:
	P_p4 > P_i4
	theta * P_p4 * nu * deltaP
	
R133:
	P_i4 > P_i4
	deltaD * P_i4 * nu * deltaP

R134:
	P_u5 > P_i5
	theta * P_u5 * nu * deltaP
    
R135:
	P_c5 > P_i5
	theta * P_c5 * nu * deltaP
	
R136:
	P_p5 > P_i5
	theta * P_p5 * nu * deltaP
	
R137:
	P_i5 > P_i5
	deltaD * P_i5 * nu * deltaP

R138:
	P_u6 > P_i6
	theta * P_u6 * nu * deltaP
    
R139:
	P_c6 > P_i6
	theta * P_c6 * nu * deltaP
	
R140:
	P_p6 > P_i6
	theta * P_p6 * nu * deltaP
	
R141:
	P_i6 > P_i6
	deltaD * P_i6 * nu * deltaP











########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1
N_u3 = 1
N_u4 = 1
N_u5 = 1
N_u6 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0
N_c6 = 0

D_u = 1
D_c = 0

P_u1 = 3
P_u2 = 3
P_u3 = 3
P_u4 = 3
P_u5 = 3
P_u6 = 3

P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0
P_c6 = 0

P_p1 = 0
P_p2 = 0
P_p3 = 0
P_p4 = 0
P_p5 = 0
P_p6 = 0

P_i1 = 0
P_i2 = 0
P_i3 = 0
P_i4 = 0
P_i5 = 0
P_i6 = 0

Acquisition = 0
Acquisition2 = 0
Deaths = 0

# Contact Rates and Contamination Probabilities #
rho_N = 3.973 # nurse direct care tasks per patient per hour
rho_D = 0.181 # doctor direct care tasks per patient per hour 
sigma = 0.054 # hand contamination probability
psi = 0.0527 # successful colonization of an uncolonized patient probability

# Exit (death/discharge) rates
theta = 0.00949 # probability of death/discharge

# Admission Proportions
nu = 0.0779 # proportion of admissions of colonized with MRSA

# Handwashing and Gown/Glove Change Rates
iota_N = 6.404 #11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 1.748 #3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 2.728 #3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.744 #0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083 # natural decolonization rate median 20 days per Star*ICU trial
epsilonE = 0.125  # effectiveness of chg
epsilonU = 1.0    # rate protected people are returned to uncolonized
epsilonB = 0.04167    # rate people are bathed in chg

deltaI = 0.000343  # rate colonized become infected
deltaD = 0.000403    # rate infected die
deltaC = 0.002976      # rate infected cured
deltaP = 0.32      # % infected from colonized



