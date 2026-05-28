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
alpha_mtd = 2.0   
Kd_mtd = 40.0     
k_elim_mtd = 0.7  
n_mtd = 2.0       

# MCT (Metronomic Chemotherapy - daily low dose)
dose_mct = 5.0    
alpha_mct = 0.6   
Kd_mct = 10.0     
k_elim_mct = 0.2  
n_mct = 1.0       

# ==========================================
# 3. DIFFERENTIAL EQUATION SOLVER (3 SCENARIOS)
# ==========================================
# Tumor Volumes
V_seq = np.zeros(len(time)); V_seq[0] = V0       # 1. Sequential Therapy
V_swi = np.zeros(len(time)); V_swi[0] = V0       # 2. Intermittent (Chemo-Switch)
V_con = np.zeros(len(time)); V_con[0] = V0       # 3. Concurrent Therapy

C_mtd_seq = np.zeros(len(time)); C_mtd_seq[0] = dose_mtd
C_mct_seq = np.zeros(len(time))

C_mtd_swi = np.zeros(len(time)); C_mtd_swi[0] = dose_mtd
C_mct_swi = np.zeros(len(time))

C_mtd_con = np.zeros(len(time)); C_mtd_con[0] = dose_mtd
C_mct_con = np.zeros(len(time)); C_mct_con[0] = dose_mct

for i in range(1, len(time)):
    t = time[i]
    t_r = round(t, 1) # Float precision handling
    
    # --- A. DOSING LOGIC ---
        
    # 1. Sequential (Only MTD for first 63 days, then only MCT)
    if t_r <= 63 and t_r % 21 == 0:
        C_mtd_seq[i-1] += dose_mtd
    elif t_r > 63 and t_r % 1 == 0:
        C_mct_seq[i-1] += dose_mct
        
    # 2. Intermittent / Chemo-Switch (MTD every 21 days, MCT on rest days)
    if t_r % 21 == 0:
        C_mtd_swi[i-1] += dose_mtd
    elif t_r % 1 == 0:
        C_mct_swi[i-1] += dose_mct
        
    # 3. Concurrent (Both drugs on their own schedule without exception)
    if t_r % 21 == 0: C_mtd_con[i-1] += dose_mtd
    if t_r % 1 == 0: C_mct_con[i-1] += dose_mct

    # --- B. PHARMACOKINETICS (Elimination) ---
    def elimination(C, k_elim): return -k_elim * C * dt
    
    
    C_mtd_seq[i] = C_mtd_seq[i-1] + elimination(C_mtd_seq[i-1], k_elim_mtd)
    C_mct_seq[i] = C_mct_seq[i-1] + elimination(C_mct_seq[i-1], k_elim_mct)
    
    C_mtd_swi[i] = C_mtd_swi[i-1] + elimination(C_mtd_swi[i-1], k_elim_mtd)
    C_mct_swi[i] = C_mct_swi[i-1] + elimination(C_mct_swi[i-1], k_elim_mct)
    
    C_mtd_con[i] = C_mtd_con[i-1] + elimination(C_mtd_con[i-1], k_elim_mtd)
    C_mct_con[i] = C_mct_con[i-1] + elimination(C_mct_con[i-1], k_elim_mct)

    # --- C. PHARMACODYNAMICS (Emax / Hill Equation) ---
    def effect(C, alpha, Kd, n): return alpha * (C**n) / (Kd**n + C**n)
    
    E_seq = effect(C_mtd_seq[i], alpha_mtd, Kd_mtd, n_mtd) + effect(C_mct_seq[i], alpha_mct, Kd_mct, n_mct)
    E_swi = effect(C_mtd_swi[i], alpha_mtd, Kd_mtd, n_mtd) + effect(C_mct_swi[i], alpha_mct, Kd_mct, n_mct)
    E_con = effect(C_mtd_con[i], alpha_mtd, Kd_mtd, n_mtd) + effect(C_mct_con[i], alpha_mct, Kd_mct, n_mct)

    # --- D. RESISTANCE LOGIC (Combinations use Immune Threshold) ---
    
    res_seq = max(0, (V_seq[i-1] - (V_immune - 10)) / max(1, V_seq[i-1]))
    res_swi = max(0, (V_swi[i-1] - (V_immune - 10)) / max(1, V_swi[i-1]))
    res_con = max(0, (V_con[i-1] - (V_immune - 10)) / max(1, V_con[i-1]))

    # --- E. TUMOR VOLUME UPDATE ---
    def dV(V, E, res): return r * V * (1 - V / K) - (E * res * V)
    
    # Combinations can reach 0 (Immune Cure)
    V_seq[i] = max(0, V_seq[i-1] + dV(V_seq[i-1], E_seq, res_seq) * dt)
    V_swi[i] = max(0, V_swi[i-1] + dV(V_swi[i-1], E_swi, res_swi) * dt)
    V_con[i] = max(0, V_con[i-1] + dV(V_con[i-1], E_con, res_con) * dt)

    # --- F. IMMUNOLOGICAL CURE CHECK ---
    if V_seq[i] <= V_immune: V_seq[i] = 0.0
    if V_swi[i] <= V_immune: V_swi[i] = 0.0
    if V_con[i] <= V_immune: V_con[i] = 0.0

# ==========================================
# 5. VISUALIZATION: COMBINATION SYNERGY (OPTIMIZED)
# ==========================================
fig2, ax3 = plt.subplots(figsize=(14, 7))

# Plot Combination Therapies
ax3.plot(time, V_seq, label="1. Sequential Therapy (63 Days MTD -> MCT Maintenance)", color="#e67e22", linewidth=2.5, linestyle="--")
ax3.plot(time, V_swi, label="2. Intermittent Therapy (Chemo-Switch)", color="#8e44ad", linewidth=2.5, linestyle="--")
ax3.plot(time, V_con, label="3. Concurrent Therapy (Simultaneous)", color="#27ae60", linewidth=3)

# Threshold Lines and Shading
ax3.axhline(y=V_dormancy, color='gray', linestyle='--', linewidth=1.5, label=f"Dormancy Threshold ({V_dormancy} mm³)")
ax3.axhline(y=V_immune, color='green', linestyle=':', linewidth=2, label=f"Immune Cure Threshold ({V_immune} mm³)")

# VISUAL ENHANCEMENT: Shade the Immune Cure Zone
ax3.fill_between(time, 0, V_immune, color='green', alpha=0.1, label="Immunological Clearance Zone")

# VISUAL ENHANCEMENT: Marker for Sequential Therapy transition
ax3.axvline(x=63, color="#225ae6", linestyle='-', alpha=0.7)
ax3.text(64, 1100, "<- Sequential: MTD Ends, MCT Starts", color="#225ae6", fontsize=12, fontstyle='italic', horizontalalignment='left', verticalalignment='top')

# Aesthetics and Labels
ax3.set_title("Tumor Dynamics: Synergy of Combination Therapies", fontsize=15, fontweight="bold")
ax3.set_xlabel("Time (Days)", fontsize=13)
ax3.set_ylabel("Tumor Volume (mm³)", fontsize=13)

# VISUAL ENHANCEMENT: Limit X-axis to 150 to zoom in on the cure events
ax3.set_xlim(0, 100)
ax3.set_ylim(0, 1200)

ax3.grid(alpha=0.4, linestyle=':')
ax3.legend(fontsize=12, loc="upper left")

plt.tight_layout()
plt.show()
