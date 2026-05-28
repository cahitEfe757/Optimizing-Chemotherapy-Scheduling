import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PATIENT AND TUMOR PARAMETERS
# ==========================================
V0 = 1000.0         # Initial tumor volume
K = 2500.0          # Maximum tumor carrying capacity
V_dormancy = 200.0  # Dormancy Threshold (Barrier for monotherapies)
V_immune = 50.0     # Immunological Cure Threshold: Immune system clears tumor below this
r = 0.08            # Intrinsic tumor growth rate
days = 180          # Total simulation days
dt = 0.1
time = np.arange(0, days, dt)

# ==========================================
# 2. DRUG PARAMETERS
# ==========================================
# MTD (Maximum Tolerated Dose - every 21 days)
dose_mtd = 100.0  
alpha_mtd = 5.0   
Kd_mtd = 40.0     
k_elim_mtd = 0.5  
n_mtd = 2.0       

# MCT (Metronomic Chemotherapy - daily low dose)
dose_mct = 5.0    
alpha_mct = 0.6   
Kd_mct = 10.0     
k_elim_mct = 0.2  
n_mct = 1.0       

# ==========================================
# 3. DIFFERENTIAL EQUATION SOLVER 
# ==========================================
# Tumor Volumes
V_mtd = np.zeros(len(time)); V_mtd[0] = V0       # Monotherapy MTD
V_mct = np.zeros(len(time)); V_mct[0] = V0       # Monotherapy MCT

# Drug Concentrations
C_mtd_alone = np.zeros(len(time)); C_mtd_alone[0] = dose_mtd
C_mct_alone = np.zeros(len(time)); C_mct_alone[0] = dose_mct


for i in range(1, len(time)):
    t = time[i]
    t_r = round(t, 1) # Float precision handling
    
    # --- A. DOSING LOGIC ---
    # Monotherapy MTD and MCT
    if t_r % 21 == 0: C_mtd_alone[i-1] += dose_mtd
    if t_r % 1 == 0: C_mct_alone[i-1] += dose_mct
        
    # --- B. PHARMACOKINETICS (Elimination) ---
    def elimination(C, k_elim): return -k_elim * C * dt
    
    C_mtd_alone[i] = C_mtd_alone[i-1] + elimination(C_mtd_alone[i-1], k_elim_mtd)
    C_mct_alone[i] = C_mct_alone[i-1] + elimination(C_mct_alone[i-1], k_elim_mct)

    # --- C. PHARMACODYNAMICS (Emax / Hill Equation) ---
    def effect(C, alpha, Kd, n): return alpha * (C**n) / (Kd**n + C**n)
    
    E_mtd_alone = effect(C_mtd_alone[i], alpha_mtd, Kd_mtd, n_mtd)
    E_mct_alone = effect(C_mct_alone[i], alpha_mct, Kd_mct, n_mct)

    # --- D. RESISTANCE LOGIC (Combinations use Immune Threshold) ---
    res_mtd = max(0, (V_mtd[i-1] - V_dormancy) / max(1, V_mtd[i-1]))
    res_mct = max(0, (V_mct[i-1] - V_dormancy) / max(1, V_mct[i-1]))
    

    # --- E. TUMOR VOLUME UPDATE ---
    def dV(V, E, res): return r * V * (1 - V / K) - (E * res * V)
    
    # Monotherapies cannot breach V_dormancy (200)
    V_mtd[i] = max(V_dormancy, V_mtd[i-1] + dV(V_mtd[i-1], E_mtd_alone, res_mtd) * dt)
    V_mct[i] = max(V_dormancy, V_mct[i-1] + dV(V_mct[i-1], E_mct_alone, res_mct) * dt)
    


# ==========================================
# 4. VISUALIZATION 1: MONOTHERAPY FAILURE
# ==========================================
fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# TOP PANEL: Drug Concentrations (PK)
ax1.plot(time, C_mtd_alone, label="MTD PK Profile (21-Day Cycle)", color="#e74c3c", linewidth=2)
ax1.plot(time, C_mct_alone, label="MCT PK Profile (Daily Continuous)", color="#3498db", linewidth=2)
ax1.set_title("Pharmacokinetic Profile: Monotherapies", fontsize=14, fontweight="bold")
ax1.set_ylabel("Drug Concentration", fontsize=12)
ax1.grid(alpha=0.4, linestyle=':')
ax1.legend(fontsize=11, loc="upper right")

# BOTTOM PANEL: Tumor Volume (PD)
ax2.plot(time, V_mtd, label="MTD Monotherapy (Relapse Cycle)", color="#c0392b", linewidth=2.5)
ax2.plot(time, V_mct, label="MCT Monotherapy (Dormancy Plateau)", color="#2980b9", linewidth=2.5)

ax2.axhline(y=V_dormancy, color='gray', linestyle='--', linewidth=1.5, label=f"Dormancy Threshold ({V_dormancy} mm³)")
ax2.axhline(y=V_immune, color='green', linestyle=':', linewidth=2.5, label=f"Immune Cure Threshold ({V_immune} mm³) - Target")

ax2.set_title("Tumor Dynamics: Monotherapy Limitations", fontsize=14, fontweight="bold")
ax2.set_xlabel("Time (Days)", fontsize=13)
ax2.set_ylabel("Tumor Volume (mm³)", fontsize=13)
ax2.set_ylim(0, 1200)
ax2.grid(alpha=0.4, linestyle=':')
ax2.legend(fontsize=11, loc="upper right")

plt.tight_layout()
plt.show()

