
# Imports 
import numpy as np
from scipy.optimize import fsolve
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle


# units 
mm = 0.001
m = 1
N = .001
kN = 1
kPa = 1

# Inputs 
H = 8000 *mm 
t = 240 *mm  
s = 1200 *mm 
db = 25 *mm
d = t/2 
fblock= 25  

# Wall Properties

faim = 0.6
fais = 0.85
emu = 0.003
k= 1

lambda_h = k * H /t 




def cross_section(t, s, db,fblock):
    
    if db == 10 * mm:
        As = 100 * mm**2
    elif db == 15 * mm:
        As = 200 * mm**2
    elif db == 20 * mm:
        As = 300 * mm**2
    elif db == 25 * mm:
        As = 500 * mm**2


    
    # Material properties

    if fblock == 10:
        fm_g = 5
        fm_ug = 6.5
    elif fblock == 15:
        fm_g = 7.5
        fm_ug = 10
    elif fblock == 20:
        fm_g = 10
        fm_ug = 13
    elif fblock == 25:
        fm_g = 11.75
        fm_ug = 15.25
    elif fblock == 30:
        fm_g = 13.5
        fm_ug = 17.5

    if t == 140 *mm: 
        tf = 26 *mm 
        rho_g = 3.03 *kN/m**2
        rho_ug = 1.67 *kN/m**2
    elif t == 190 *mm:  
        tf = 36.2 *mm
        rho_g = 4.12 *kN/m**2
        rho_ug = 2.19 *kN/m**2
    elif t == 240 *mm:  
        tf = 38.6 *mm
        rho_g = 5.22 *kN/m**2
        rho_ug = 2.62 *kN/m**2
    elif t == 290 *mm:  
        tf = 41.4 *mm
        rho_g = 6.32 *kN/m**2
        rho_ug = 3.059 *kN/m**2
    
    beff = np.minimum(4*t, s)
    beff_m_1= 1000 *mm
    beff_m_2= beff * 1/s
    Aseff_m = As *1/s
    bg_m = 200*mm * 1/s
    bug_m_1 = 1000*mm - bg_m  # Bars in compression 
    bug_m_2 = beff_m_2- bg_m  #Bars in tension

    # Area effective 
    A_gr= bg_m * t
    A_ug_1 = bug_m_1 * 2*tf    # Bars in compression 
    A_ug_2 = bug_m_2 * 2*tf   # Bars in tension

    Ae_1 = A_gr + A_ug_1 # Bars in compression 
    Ae_2 = A_gr + A_ug_2 # Bars in tension

    # fm effective

    fm_e_1 = (fm_g * A_gr + fm_ug * A_ug_1) / (Ae_1)
    fm_e_2 = (fm_g * A_gr + fm_ug * A_ug_2) / (Ae_2)
    E_m = 850 * fm_e_2

    # Inertia 
    I_gross_gr=  beff_m_1 * t**3 / 12
    I_gross_ug_1 = (beff_m_1 * tf**3 / 12 + beff_m_1 * tf * (t/2 - tf/2)**2) *2
    I_gross_eff = (I_gross_gr * A_gr + + I_gross_ug_1 * A_ug_1)/ (Ae_1)  
    S_eff = I_gross_eff / (t/2)
    n = 200000/E_m
    kd = sp.Symbol('kd', real=True)
    eq = n * Aseff_m * (d - kd) - (beff_m_2 * kd**2) / 2
    solutions  = sp.solve(eq, kd)
    kd= [sol.evalf() for sol in solutions if sol.is_real and sol > 0][0] 
    I_cr_eff = n * Aseff_m * (d- kd)**2 + (beff_m_2 * kd**3) / 3    
    ek = S_eff / Ae_1

   # Self Weight 
    rho_SW = (rho_g * A_gr + rho_ug * A_ug_1) / (Ae_1) 

    return beff_m_1, beff_m_2, Aseff_m, bg_m, bug_m_1, bug_m_2, Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf 

beff_m_1, beff_m_2, Aseff_m, bg_m, bug_m_1, bug_m_2, Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n ,E_m , ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf =cross_section(t, s,db,fblock)
    

# Applied Loads 
P_DL = 10 *kN
P_LL = 0 *kN 
e_P = 0 *mm
W = 1 * kPa

e_P = np.maximum(0.1 * t, e_P)
# Self Weight
P_SW_mid = rho_SW  * H/2  # at mid span 
P_SW_base = rho_SW  * H # at base

# moment due to wind 
M_lateral = W * H**2 / 8
M_unf = (P_DL + P_LL) * e_P + M_lateral

M_lateral_F1 = 0
M_lateral_F2 = 0
M_lateral_F3 = 0
M_lateral_F4 = 1.4 *  M_lateral
M_lateral_F5 = 1.4 *  M_lateral

M_lateral_F = [M_lateral_F1, M_lateral_F2, M_lateral_F3, M_lateral_F4, M_lateral_F5]    

# Factored Primary Moment at midspan 
M_F1 = np.maximum(1.4 * P_DL * e_P /2, 1.4 * (P_DL  + P_SW_mid) * 0.1 *t)
M_F2 = np.maximum(1.25 * P_DL * e_P /2+ 1.5 * P_LL * e_P /2, 1.25 * (P_DL  + P_SW_mid) * 0.1 *t + 1.5 * P_LL * 0.1 *t)
M_F3 = np.maximum(0.9 * P_DL * e_P/2 + 1.5 * P_LL * e_P  /2, 0.9 * (P_DL  + P_SW_mid) * 0.1 *t + 1.5 * P_LL * 0.1 *t)
M_F4 = np.maximum(1.25 * P_DL * e_P /2+ 1.4 * M_lateral, 1.25 * (P_DL  + P_SW_mid) * 0.1 *t)
M_F5 = np.maximum(0.9 * P_DL * e_P /2+ 1.4 * M_lateral, 0.9 * (P_DL  + P_SW_mid) * 0.1 *t)

M_F = [M_F1, M_F2, M_F3, M_F4, M_F5]

# Factored Primary Axial Load @ midspan
P_F1 = 1.4 * (P_DL  + P_SW_mid)
P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
P_F4 = 1.25 * (P_DL  + P_SW_mid) 
P_F5 = 0.9 * (P_DL  + P_SW_mid) 

P_F = [P_F1, P_F2, P_F3, P_F4, P_F5]
# top Moment 
M1_F1 =  np.maximum(1.4 * P_DL * e_P, 1.4 * P_DL *0.1*t)
M1_F2 =  np.maximum(1.25 * P_DL * e_P+ 1.5 * P_LL * e_P, 1.25 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t) 
M1_F3 =  np.maximum(0.9 * P_DL * e_P + 1.5 * P_LL * e_P, 0.9 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t)  
M1_F4 =  np.maximum(1.25 * P_DL * e_P , 1.25 * P_DL * 0.1 *t)
M1_F5 =  np.maximum(0.9 * P_DL * e_P  , 0.9 * P_DL * 0.1 *t)

Mtop_F= [M1_F1, M1_F2, M1_F3, M1_F4, M1_F5]
# Bottom Moment 
M2_F1 = np.maximum(0,  1.4 * (P_DL  + P_SW_base) *0.1*t)
M2_F2 = np.maximum(0,  1.25 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_LL * 0.1 *t) 
M2_F3 = np.maximum(0,  0.9 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_LL * 0.1 *t)  
M2_F4 = np.maximum(0, 1.25 * (P_DL  + P_SW_base) * 0.1 *t)
M2_F5 = np.maximum(0, 0.9 * (P_DL  + P_SW_base) * 0.1 *t)

Mbot_F= [M2_F1, M2_F2, M2_F3, M2_F4, M2_F5]

# Moment Ratio
M_ratio =  np.minimum(np.array(Mtop_F), np.array(Mbot_F)) / np.maximum(np.array(Mtop_F), np.array(Mbot_F))

ev = np.array([M_F1, M_F2, M_F3, M_F4, M_F5]) / np.array([P_F1, P_F2, P_F3, P_F4, P_F5])
# Rigidty Coefficient 
betad = P_DL * e_P / M_unf
betad=0
faiE = 0.75
EI_eff_raw = E_m * (0.25 * I_gross_eff - (0.25 * I_gross_eff - I_cr_eff)* ((ev - ek) / (2 * ek))) 
EI_eff = np.clip(EI_eff_raw, E_m * I_cr_eff, 0.25 * E_m * I_gross_eff)
Rigidty_c = faiE * EI_eff / (1+0.5*betad)

Pcr = np.pi**2 * Rigidty_c/ (k*H)**2*1000
# Lateral Force Coefficient
if lambda_h < 30:
    Cm = np.where(np.array(M_lateral_F) / np.array(M_F) > 0.5,
                1,
                np.maximum(0.6 + 0.4 * M_ratio, 0.4))
else: 
    Cm = [1] * len(M_F)

MagFactor=  np.array(Cm) / (1-(np.array(P_F)/np.array(Pcr)))

# Design Moment

Mt_F = np.array(M_F) * np.array(MagFactor)



## Interaction Diagram 

# Point 1 
# Point 1 : Prmax 

def solve_betaC(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t):
    # Compute Prmax
    Prmax = 0.8 * faim * 0.85 * fm_e_1 * Ae_1 * 1000  # kN
    
    # Define the equation to solve for betaC
    def equation(betaC):
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
        F_g_bottom = faim * 0.85 * fm_ug * (betaC - (t - tf)) * bug_m_1 * 1000 if betaC >= (t - tf) else 0
        
        if betaC >= (t - tf):
            Pr = Fg + F_ug_top + F_g_bottom
        else:
            Pr = Fg + F_ug_top
        return Pr - Prmax
    
    # Solve for betaC
    betaC_solution = fsolve(equation, t)  # Initial guess is t
    betaC = betaC_solution[0]
    
    # Recalculate forces with the obtained betaC

    Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
    F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
    F_g_bottom = faim * 0.85 * fm_ug * (betaC - (t - tf)) * bug_m_1 * 1000 if betaC >= (t - tf) else 0
    
    if betaC >= (t - tf):
        Mr = Fg *(t/2 - betaC/2) + F_ug_top * (t/2-tf/2) - F_g_bottom * (t/2 - (tf-(t-betaC))/2) # moment at center 
        Mr2 = Fg *betaC/2 + F_ug_top * tf/2 + F_g_bottom * (t+(tf-(t-betaC)/2)-tf) - Prmax * t/2 # moment at top 
    else : 
        Mr = Fg * (t/2 - betaC/2) + F_ug_top * (t/2 - tf/2) # moment at center
        Mr2 = Fg *betaC/2 + F_ug_top * tf/2  - Prmax * t/2 # moment at top 
    return betaC, Fg, F_ug_top, F_g_bottom, Prmax, Mr
PMax= solve_betaC(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t)
betaC1=PMax[0]

# Steel in compression 
# def Point2(betaC1,faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t, d, num_points=3):
#     betaC_values = np.linspace(0.8*(t - d), betaC1, num=num_points)[::-1]
#     results = []
    
#     for betaC in betaC_values:
#         Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
#         F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
#         Pr = Fg + F_ug_top
#         Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2)
#         e = Mr / Pr 
#         if e > t/3 : 
#             Mr = Pr * t/3
#         results.append((betaC, Fg, F_ug_top, Pr, Mr))
    
#     return results

# Point2=Point2(betaC1,faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t, d, num_points=15)

#  # Steel in tension 
# def Point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t, d, Aseff_m, num_points):
#     betaC_values = np.linspace(tf / 2, 0.8 * (t - d), num=num_points)[::-1]
#     results = []
    
#     Mr_y, Pr_y = None, None  # Initialize to None
#     i = 0  # Initialize index
    
#     while i < len(betaC_values):
#         betaC = betaC_values[i]
#         c = betaC / 0.8
#         es = emu * (d - c) / c     
#         fs = np.minimum(200000 * es, 400)
#         Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        
#         if betaC < tf:
#             F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
#         else: 
#             F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
        
#         Ts = fais * Aseff_m * fs * 1000
#         Pr = Fg + F_ug_top - Ts 
        
#         if betaC < tf:
#             Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
#         else: 
#             Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
        
#         e = Mr / Pr 

#         # Store the first values where es >= 0.002
#         if es >= 0.002 and Mr_y is None and Pr_y is None:
#             ey, Mr_y, Pr_y = es, Mr, Pr
        
#         if Pr <= 0:
#             break  # Exit loop if Pr is non-positive

#         results.append((betaC, Fg, F_ug_top, Pr, Mr))
#         i += 1  # Increment index
    
#     return results, Mr_y, Pr_y, ey

# # Example function call
# Point3, Mr_y, Pr_y, ey = Point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t, d, Aseff_m, num_points=40)

def calculate_point2(betaC1, faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t, d, num_points=3):
    betaC_values = np.linspace(0.8*(t - d), betaC1, num=num_points)[::-1]
    results = []
    
    for betaC in betaC_values:
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
        Pr = Fg + F_ug_top
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2)
        e = Mr / Pr 
        if e > t/3:
            Mr = Pr * t/3
        results.append((betaC, Fg, F_ug_top, Pr, Mr))
    
    return results

def calculate_point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t, d, Aseff_m, num_points):
    betaC_values = np.linspace(tf / 2, 0.8 * (t - d), num=num_points)[::-1]
    results = []
    
    Mr_y, Pr_y = None, None  # Initialize to None
    i = 0  # Initialize index
    
    while i < len(betaC_values):
        betaC = betaC_values[i]
        c = betaC / 0.8
        es = emu * (d - c) / c     
        fs = np.minimum(200000 * es, 400)
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        
        if betaC < tf:
            F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
        else: 
            F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
        
        Ts = fais * Aseff_m * fs * 1000
        Pr = Fg + F_ug_top - Ts 
        
        if betaC < tf:
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
        else: 
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
        
        e = Mr / Pr 

        # Store the first values where es >= 0.002
        if es >= 0.002 and Mr_y is None and Pr_y is None:
            ey, Mr_y, Pr_y = es, Mr, Pr
        
        if Pr <= 0:
            ey=0
            break  # Exit loop if Pr is non-positive

        results.append((betaC, Fg, F_ug_top, Pr, Mr))
        i += 1  # Increment index
    
    return results, Mr_y, Pr_y, ey

# Point 4 (pure tension moment)
def calculate_pure_moment(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t):
    # Define the equation to solve for betaC
    def equation(betaC):
        Ptarget = 0.001 
        c = betaC / 0.8
        es = emu * (d - c) / c     
        fs = np.minimum(200000 * es, 400)
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        
        if betaC < tf:
            F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
        else: 
            F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
        
        Ts = fais * Aseff_m * fs * 1000
        Pr = Fg + F_ug_top - Ts 
        if betaC < tf:
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
        else:
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
        return Pr - Ptarget
    
    # Solve for betaC
    betaC_solution = fsolve(equation, tf/2)  # Initial guess is t
    betaC = betaC_solution[0]
    
    c = betaC / 0.8
    es = emu * (d - c) / c     
    fs = np.minimum(200000 * es, 400)
    Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
    
    if betaC < tf:
        F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
    else: 
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
    
    Ts = fais * Aseff_m * fs * 1000
    Pr = Fg + F_ug_top - Ts 
    if betaC < tf:
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
    else: 
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
    return betaC, Fg, F_ug_top, Pr, Mr
# def PureMoment(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t):
    
    # Define the equation to solve for betaC
    def equation(betaC):
        Ptarget = 0.001 
        c = betaC / 0.8
        es = emu * (d - c) / c     
        fs = np.minimum(200000 * es, 400)
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        
        if betaC < tf:
            F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
        else: 
            F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
        
        Ts = fais * Aseff_m * fs * 1000
        Pr = Fg + F_ug_top - Ts 
        if betaC<tf: 
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
        else:
             Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
        return Pr - Ptarget
    
    # Solve for betaC
    betaC_solution = fsolve(equation, tf/2)  # Initial guess is t
    betaC = betaC_solution[0]
    
    c = betaC / 0.8
    es = emu * (d - c) / c     
    fs = np.minimum(200000 * es, 400)
    Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
    
    if betaC < tf:
        F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
    else: 
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
    
    Ts = fais * Aseff_m * fs * 1000
    Pr = Fg + F_ug_top - Ts 
    if betaC < tf:
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
    else: 
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
    return betaC, Fg, F_ug_top, Pr, Mr
# PureMoment=PureMoment(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t)

# P= [PMax[-2]] + [PMax[-2]] + [result[3] for result in Point2] + [result[3] for result in Point3] +[ PureMoment[-2]]
# M= [0]+ [PMax[-1]] + [result[4] for result in Point2]+ [result[4] for result in Point3] + [PureMoment[-1]]

# # Create the plot
# plt.figure(figsize=(8,6))
# plt.plot(M, P,  linestyle='-', color='b', label='Moment vs. Axial Force')
# plt.scatter(Mt_F, P_F, marker='x', color='r', label='Applied Loads')

# # Labels and title
# plt.xlabel("Moment (M)")
# plt.ylabel("Axial Force (P)")
# plt.title("Moment vs. Axial Force Relationship")

# # Grid and legend
# plt.grid(True, linestyle='--', alpha=0.6)
# plt.legend()
# plt.ylim(0, 1.1 * max(P))
# plt.xlim(0, 1.1 * max(max(M), max(Mt_F)))


# Show the plot
# plt.show()



## Cros section Plot 

# def draw_blocks(height, faceshell, width, num_blocks, gap, s):
#     fig, ax = plt.subplots(figsize=(num_blocks * 4, 2))

#     # Determine the grouting pattern based on spacing s
#     grout_pattern = s // 200  # Number of cells in a grouting cycle

#     for i in range(num_blocks):
#         x_offset = i * (width + gap)

#         # Outer rectangle (block boundary)
#         outer_rect = patches.Rectangle((x_offset, 0), width, height, linewidth=2, edgecolor='black', facecolor='lightgray')
#         ax.add_patch(outer_rect)

#         # Openings (inner rectangles)
#         opening1_color = 'white'
#         opening2_color = 'white'
#         dot_positions = []

#         # Apply grouting pattern
#         if (i * 2) % grout_pattern == 0:
#             opening1_color = 'darkgrey'  # Grout the first opening
#             dot_positions.append((x_offset + width * 0.1 + width * 0.35 / 2, faceshell + (height - 2 * faceshell) / 2))
#         if (i * 2 + 1) % grout_pattern == 0:
#             opening2_color = 'darkgrey'  # Grout the second opening
#             dot_positions.append((x_offset + width * 0.55 + width * 0.35 / 2, faceshell + (height - 2 * faceshell) / 2))

#         # Draw openings
#         opening1 = patches.Rectangle((x_offset + width * 0.1, faceshell), width * 0.35, height - 2 * faceshell, linewidth=2, edgecolor='black', facecolor=opening1_color)
#         opening2 = patches.Rectangle((x_offset + width * 0.55, faceshell), width * 0.35, height - 2 * faceshell, linewidth=2, edgecolor='black', facecolor=opening2_color)

#         ax.add_patch(opening1)
#         ax.add_patch(opening2)
            
#         # Add black dots in the center of grouted cells
#         for (x, y) in dot_positions:
#             ax.plot(x, y, 'ko')
#     # Add thickness dimension text
#     ax.text(x_offset - width *num_blocks/ 3, -10, f'Thickness: {height} mm', ha='center', va='top', fontsize=15)
#     # Formatting
#     ax.set_xlim(-1, num_blocks * (width + gap))
#     ax.set_ylim(-1, height + 1)
#     ax.set_aspect('equal', adjustable='box')
#     plt.axis('off')  # Hide axes
    # plt.show()

# Example function call
# draw_blocks(height=t/mm, faceshell=tf*1.25/mm, width=400, num_blocks=5, gap=10, s=s/mm)


def draw_blocks(height, faceshell, width, num_blocks, gap, s):
    # Use a non-GUI backend
    plt.switch_backend('Agg')
    
    fig, ax = plt.subplots(figsize=(num_blocks * 4, 2))

    # Determine the grouting pattern based on spacing s
    grout_pattern = s // 200  # Number of cells in a grouting cycle

    for i in range(num_blocks):
        x_offset = i * (width + gap)

        # Outer rectangle (block boundary)
        outer_rect = Rectangle((x_offset, 0), width, height, linewidth=2, edgecolor='black', facecolor='lightgray')
        ax.add_patch(outer_rect)

        # Openings (inner rectangles)
        opening1_color = 'white'
        opening2_color = 'white'
        dot_positions = []

        # Apply grouting pattern
        if (i * 2) % grout_pattern == 0:
            opening1_color = 'darkgrey'  # Grout the first opening
            dot_positions.append((x_offset + width * 0.1 + width * 0.35 / 2, faceshell + (height - 2 * faceshell) / 2))
        if (i * 2 + 1) % grout_pattern == 0:
            opening2_color = 'darkgrey'  # Grout the second opening
            dot_positions.append((x_offset + width * 0.55 + width * 0.35 / 2, faceshell + (height - 2 * faceshell) / 2))

        # Draw openings
        opening1 = Rectangle((x_offset + width * 0.1, faceshell), width * 0.35, height - 2 * faceshell, linewidth=2, edgecolor='black', facecolor=opening1_color)
        opening2 = Rectangle((x_offset + width * 0.55, faceshell), width * 0.35, height - 2 * faceshell, linewidth=2, edgecolor='black', facecolor=opening2_color)

        ax.add_patch(opening1)
        ax.add_patch(opening2)
            
        # Add black dots in the center of grouted cells
        for (x, y) in dot_positions:
            ax.plot(x, y, 'ko')
    
    # Add thickness dimension text
    ax.text(-width * num_blocks / 3, -10, f'Thickness: {height} mm', ha='center', va='top', fontsize=15)
    
    # Formatting
    ax.set_xlim(-1, num_blocks * (width + gap))
    ax.set_ylim(-1, height + 1)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')  # Hide axes

    # Return the figure instead of displaying it
    return fig
