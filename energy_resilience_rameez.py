import numpy as np
from scipy.optimize import linprog
from scipy.io import savemat # ADDED: Import for saving to MATLAB format
import pandas as pd
import matplotlib.pyplot as plt 
from typing import Tuple

T = 672 # February 2025 (28 days * 24 hours)

# Economic Parameters
C_Thermal = 50 # Thermal Plant Marginal Cost (£/MWh)
C_VoLL = 5000 # (Assumed 100 times C_Thermal) Value of Lost Load (VoLL) / Cost of unserved energy (£/MWh)

# System Component Capacities (MW / MWh)
P_Th_Max = 250 #Thermal max capacity
P_Wind_Max = 150 
Line_Max = 100
P_Batt_Rate = 100
E_Batt_Cap = 2400 # 100 MW * 24h capacity

# Battery Parameters (State of Charge (SOC))
Eff_Ch = 0.95 # Charging Efficiency (85% to 95% for Sulphur-flow battery)
Eff_Dis = 0.95 # Discharging Efficiency (85% to 95% for Sulphur-flow battery)
E_Batt_Init = E_Batt_Cap * 0.5 # Initial State of Charge (50%)
E_Batt_Min = E_Batt_Cap * 0.1 # Minimum State of Charge (10%)

# --- Data Loading (Load and Wind Profiles) ---
data = pd.read_csv('Hourly Profiles.csv')
# Scale load profile by the 100 MW peak demand
P_Load_Profile = data.iloc[0:T, 2].values * 100 #Used directly as parameter
P_Wind_Profile = data.iloc[0:T, 1].values #Used as a multiplier fraction.


# === 2. LINEAR PROGRAMMING SETUP ===

# Decision Variable Vector X (12 * T variables):
# [P_Th; P_Curtail; P_Dis; P_Ch; P_Unserved; F_12+; F_12-; F_13+; F_13-; F_23+; F_23-; E_SOC]
N_vars = 12 * T #Over each hour of total time horizon

# Objective Function (Minimize Total Cost)
C = np.zeros(N_vars)
C[0:T] = C_Thermal 	 # Cost of Thermal Power (P_Th)
C[4 * T: 5 * T] = C_VoLL # Cost of Unserved Load (VoLL penalty) (P_unserved)
# All other variables have a default cost of zero.

# === 3. EQUALITY CONSTRAINTS (Aeq * X = Beq) ===
#Aeq : Coefficient matrix
#Beq : Constants vector
#X   : Decision variable

#Kirchoff's current rule
# Nodal Balance (3*T) + SOC Balance (T) = 4*T rows
N_eq = 4 * T
Aeq = np.zeros((N_eq, N_vars))
Beq = np.zeros(N_eq)
row = 0 #counter variable

# Line Flow Block Start Indices (for reference): Start for P_Th is 0T to 1T
# F_12+: 5T | F_12-: 6T | F_13+: 7T | F_13-: 8T | F_23+: 9T | F_23-: 10T | E_SOC: 11T

# 3.1 Nodal Power Balance (3*T rows) 
# The transport model is used, where F_ij = F_ij+ - F_ij-. #Accounts for bidirectionality
for t in range(T):
	t_idx = t

	# Bus 1 (Thermal Generator): P_Th - F_12 - F_13 = 0 (Bus Lines are outputs)
	# P_Th - (F_12+ - F_12-) - (F_13+ - F_13-) = 0
	Aeq[row, t_idx] = 1
	Aeq[row, 5 * T + t_idx] = -1; Aeq[row, 6 * T + t_idx] = 1 	# -F_12+ + F_12-
	Aeq[row, 7 * T + t_idx] = -1; Aeq[row, 8 * T + t_idx] = 1 	# -F_13+ + F_13-
	Beq[row] = 0
	row += 1

	# Bus 2 (Wind + Battery): P_Dis - P_Ch + F_12 - F_23 + P_Curtail = P_Wind_Avail
	# P_Dis - P_Ch + P_Curtail + (F_12+ - F_12-) - (F_23+ - F_23-) = P_Wind_Avail
	Aeq[row, T + t_idx] = 1 	 # P_Curtail (Block 2)
	Aeq[row, 2 * T + t_idx] = 1 	 # P_Discharge (Block 3)
	Aeq[row, 3 * T + t_idx] = -1 	 # P_Charge (Block 4)
	
	Aeq[row, 5 * T + t_idx] = 1; Aeq[row, 6 * T + t_idx] = -1 	 # +F_12+ - F_12- (Flow 1->2 is positive)
	Aeq[row, 9 * T + t_idx] = -1; Aeq[row, 10 * T + t_idx] = 1 	# -F_23+ + F_23- (Flow 2->3 is positive)
	Beq[row] = P_Wind_Profile[t_idx] * P_Wind_Max # P_Wind_Avail (RHS) Line 31 (Multiplier)
	row += 1

	# Bus 3 (Load): P_Unserved - F_13 - F_23 = P_Load
	# P_Unserved - (F_13+ - F_13-) - (F_23+ - F_23-) = P_Load
	Aeq[row, 4 * T + t_idx] = 1 	 # P_Unserved (Block 5)
	
	Aeq[row, 7 * T + t_idx] = -1; Aeq[row, 8 * T + t_idx] = 1 	# -F_13+ + F_13- (Flow 1->3 into Bus 3)
	Aeq[row, 9 * T + t_idx] = -1; Aeq[row, 10 * T + t_idx] = 1 	# -F_23+ + F_23- (Flow 2->3 into Bus 3)
	Beq[row] = P_Load_Profile[t_idx] # P_Load (RHS) Line 30
	row += 1

# 3.2 Battery SOC Balance (T rows)
# E_SOC(t) - E_SOC(t-1) - Eff_Ch*P_Ch(t) + (1/Eff_Dis)*P_Dis(t) = 0
for t in range(T):
	t_idx = t

	Aeq[row, 11 * T + t_idx] = 1 # +E_SOC(t)

	if t > 0:
		Aeq[row, 11 * T + t_idx - 1] = -1 # -E_SOC(t-1)
	else:
		# Initial condition E_SOC(0) = E_Batt_Init
		Beq[row] += E_Batt_Init

	Aeq[row, 3 * T + t_idx] = -Eff_Ch 	 # -Eff_Ch * P_Ch(t)
	Aeq[row, 2 * T + t_idx] = 1 / Eff_Dis # +(1/Eff_Dis) * P_Dis(t)
	row += 1

# === 4. INEQUALITY CONSTRAINTS (A * X <= B) ===

# Thermal Max (T) + SOC Min/Max (2*T) + Batt Rate Limits (2*T) = 5*T rows
N_ineq = 5 * T
A = np.zeros((N_ineq, N_vars))
B = np.zeros(N_ineq)
row = 0

# 4.1 Thermal Generation Limit (P_Th <= P_Th_Max)
for t in range(T):
	A[row, t] = 1
	B[row] = P_Th_Max
	row += 1

# 4.2 Battery SOC Max (E_SOC <= E_Batt_Cap)
for t in range(T):
	A[row, 11 * T + t] = 1
	B[row] = E_Batt_Cap
	row += 1

# 4.3 Battery SOC Min (E_SOC >= E_Batt_Min) => -E_SOC <= -E_Batt_Min
for t in range(T):
	A[row, 11 * T + t] = -1
	B[row] = -E_Batt_Min
	row += 1

# 4.4 Battery Discharge Rate (P_Dis <= P_Batt_Rate)
for t in range(T):
	A[row, 2 * T + t] = 1
	B[row] = P_Batt_Rate
	row += 1

# 4.5 Battery Charge Rate (P_Ch <= P_Batt_Rate)
for t in range(T):
	A[row, 3 * T + t] = 1
	B[row] = P_Batt_Rate
	row += 1

# === 5. VARIABLE BOUNDS ===
Lb = np.zeros(N_vars)
Ub = np.ones(N_vars) * np.inf

# Line Flow Limits (F_ij+ and F_ij- <= Line_Max)
# These apply to blocks 5T through 10T, covering F_12+, F_12-, F_13+, F_13-, F_23+, F_23-
flow_blocks_start = 5*T
flow_blocks_end = 11*T
Ub[flow_blocks_start: flow_blocks_end] = Line_Max


# === 6. SCENARIO SOLVER FUNCTION ===
def run_scenario(scenario_name: str, P_Load_in: np.ndarray, P_Wind_Avail_in: np.ndarray, P_Th_Max_in: np.ndarray, A_base: np.ndarray, B_base: np.ndarray, Aeq_base: np.ndarray, Beq_base: np.ndarray, C: np.ndarray, Lb: np.ndarray, Ub: np.ndarray, T: int) -> Tuple[float, np.ndarray, np.ndarray]:
	"""
	Runs the LP model for a given scenario by updating the right-hand side (RHS)
	of the constraints (B and Beq) and solving the problem.
	"""
	scenario_B = np.copy(B_base)
	scenario_Beq = np.copy(Beq_base)

	# 1. Update Thermal Max Capacity (Inequality B vector)
	# Rows 0 to T-1 are P_Th(t) <= P_Th_Max
	scenario_B[0:T] = P_Th_Max_in

	# 2. Update Load Demand (Bus 3 Balance - Beq)
	for t in range(T):
		# Bus 3 is the 3rd row in every 3-row nodal balance block: index 3*t + 2
		scenario_Beq[3 * t + 2] = P_Load_in[t]

	# 3. Update Wind Availability (Bus 2 Balance - Beq)
	for t in range(T):
		# Bus 2 is the 2nd row in every 3-row nodal balance block: index 3*t + 1
		scenario_Beq[3 * t + 1] = P_Wind_Avail_in[t]

	# --- Solve LP ---
	bounds = list(zip(Lb, Ub))

	try:
		# Use HiGHS solver 
		res = linprog(C, A_ub=A_base, b_ub=scenario_B, A_eq=Aeq_base, b_eq=scenario_Beq, bounds=bounds, method='highs')
	except Exception:
		# Fallback to a slower method 'revised simplex' if 'highs' is not available
		res = linprog(C, A_ub=A_base, b_ub=scenario_B, A_eq=Aeq_base, b_eq=scenario_Beq, bounds=bounds, method='revised simplex')

	if res.success:
		X = res.x
		fval = res.fun
		# P_Unserved is Block 5 (indices 4*T to 5*T)
		P_Unserved = X[4 * T: 5 * T]

		print(f"\n--- Scenario: {scenario_name} ---")
		print(f"Status: Success (Cost: £{fval:,.2f} | Total Unserved: {np.sum(P_Unserved):,.2f} MWh)")
		return fval, P_Unserved, X # Also return X for potential future analysis
	else:
		print(f"\n--- Scenario: {scenario_name} ---")
		print(f"Status: Failed - {res.message}")
		return np.nan, np.zeros(T), np.zeros(N_vars)


# === 7. SCENARIO DEFINITION AND EXECUTION ===

# Dunkelflaute period (Week 2: Hour 168 to 335, inclusive)
h_start_idx = 7 * 24
h_end_idx = 14 * 24 # Slice end is exclusive (index 336 for hour 336)

P_Load_Base = P_Load_Profile
P_Wind_Avail_Base = P_Wind_Profile * P_Wind_Max
P_Th_Max_Vector = np.ones(T) * P_Th_Max

# --- SCENARIO 1: NORMAL CONDITIONS ---
C1, PU1, X1 = run_scenario('1. Normal Conditions', P_Load_Base, P_Wind_Avail_Base, P_Th_Max_Vector, A, B, Aeq, Beq, C, Lb, Ub, T)

# --- SCENARIO 2: ABNORMAL DUNKELFLAUTE ---
P_Load_Dunkelflaute = np.copy(P_Load_Base)
P_Wind_Dunkelflaute = np.copy(P_Wind_Avail_Base)

# Apply 40% load increase and 60% wind reduction during Week 2
P_Load_Dunkelflaute[h_start_idx:h_end_idx] *= 1.40
P_Wind_Dunkelflaute[h_start_idx:h_end_idx] *= (1 - 0.60) # 40% availability

C2, PU2, X2 = run_scenario('2. Abnormal Dunkelflaute', P_Load_Dunkelflaute, P_Wind_Dunkelflaute, P_Th_Max_Vector, A, B, Aeq, Beq, C, Lb, Ub, T)

# --- SCENARIO 3 & 4: THERMAL OUTAGE DURING DUNKELFLAUTE ---
P_Th_Max_outage = np.copy(P_Th_Max_Vector)
# Outage window: 48h into Week 2, lasting 24h
h_outage_start_idx = h_start_idx + 48 # Index 168 + 48 = 216
h_outage_end_idx = h_outage_start_idx + 24 # Index 240
P_Th_Max_outage[h_outage_start_idx:h_outage_end_idx] = 0 # Thermal capacity set to zero

# Assumption: In this deterministic LP model, the solution for planned (S3) and unplanned (S4)
# is identical as the constraints are the same.
C3, PU3, X3 = run_scenario('3. Planned 24h Thermal Outage', P_Load_Dunkelflaute, P_Wind_Dunkelflaute, P_Th_Max_outage, A, B, Aeq, Beq, C, Lb, Ub, T)
C4, PU4, X4 = run_scenario('4. Unplanned 24h Thermal Outage', P_Load_Dunkelflaute, P_Wind_Dunkelflaute, P_Th_Max_outage, A, B, Aeq, Beq, C, Lb, Ub, T)

# === 8. MATLAB DATA EXPORT ===
# Prepare data dictionary for MATLAB (.mat) format.
# All scenario inputs and the full decision variable vectors (X) are saved.
mat_data = {
	# Time horizon and indices
	'T': T,
	'h_start_idx': h_start_idx,
	'h_end_idx': h_end_idx,
	'h_outage_start_idx': h_outage_start_idx,
	'h_outage_end_idx': h_outage_end_idx,

	# Input Profiles and Capacity Vectors
	'P_Load_Base': P_Load_Base,
	'P_Wind_Avail_Base': P_Wind_Avail_Base,
	'P_Load_Dunkelflaute': P_Load_Dunkelflaute,
	'P_Wind_Dunkelflaute': P_Wind_Dunkelflaute,
	'P_Th_Max_Vector': P_Th_Max_Vector,
	'P_Th_Max_outage': P_Th_Max_outage,

	# SCENARIO 1: Normal Conditions
	'X1': X1,
	'Cost1': C1,
	'P_Unserved1': PU1,

	# SCENARIO 2: Abnormal Dunkelflaute
	'X2': X2,
	'Cost2': C2,
	'P_Unserved2': PU2,

	# SCENARIO 3: Planned Outage
	'X3': X3,
	'Cost3': C3,
	'P_Unserved3': PU3,

	# SCENARIO 4: Unplanned Outage
	'X4': X4,
	'Cost4': C4,
	'P_Unserved4': PU4,
	
	# Decision Variable Structure Description
	# [P_Th; P_Curtail; P_Dis; P_Ch; P_Unserved; F_12+; F_12-; F_13+; F_13-; F_23+; F_23-; E_SOC]
	'N_vars_per_T': 12,
	'Var_Index_P_Th': [0*T, 1*T-1],
	'Var_Index_P_Curtail': [1*T, 2*T-1],
	'Var_Index_P_Dis': [2*T, 3*T-1],
	'Var_Index_P_Ch': [3*T, 4*T-1],
	'Var_Index_P_Unserved': [4*T, 5*T-1],
	'Var_Index_F_12_P': [5*T, 6*T-1],
	'Var_Index_F_12_N': [6*T, 7*T-1],
	'Var_Index_F_13_P': [7*T, 8*T-1],
	'Var_Index_F_13_N': [8*T, 9*T-1],
	'Var_Index_F_23_P': [9*T, 10*T-1],
	'Var_Index_F_23_N': [10*T, 11*T-1],
	'Var_Index_E_SOC': [11*T, 12*T-1],
}

# Save the dictionary to a .mat file
savemat('LP_Optimization_Results.mat', mat_data)
print("\n--- Data Export Complete ---")
print("Results saved to 'LP_Optimization_Results.mat' for MATLAB analysis.")


# === 9. VISUALIZATION AND FINAL METRICS SUMMARY (Original Section 8) ===

'''
# === 8. VISUALIZATION AND FINAL METRICS SUMMARY ===

# 8.1 Data Preparation for Plotting
Time = np.arange(1, T + 1)
Time_Week2 = Time[h_start_idx:h_end_idx]

# 8.2 Plotting Results
plt.figure(figsize=(14, 10))

# Subplot 1: Total Unserved Energy Across All Scenarios (Full Horizon)
plt.subplot(2, 1, 1)
plt.plot(Time, PU1, label='1. Normal Conditions')
plt.plot(Time, PU2, label='2. Dunkelflaute Baseline')
plt.plot(Time, PU3, '--', label='3. Planned Outage')
plt.plot(Time, PU4, ':', label='4. Unplanned Outage')
plt.title('Unserved Energy (Resilience Metric) Across 28 Days (672 Hours)')
plt.ylabel('Unserved Energy (MW)')
plt.xlabel('Hour (t)')
plt.legend(loc='upper left')
plt.grid(True, linestyle='--')

# Subplot 2: Detailed View During Dunkelflaute Week (Week 2)
plt.subplot(2, 1, 2)
plt.plot(Time_Week2, PU2[h_start_idx:h_end_idx], 'k-', linewidth=2, label='2. Dunkelflaute Baseline')
plt.plot(Time_Week2, PU3[h_start_idx:h_end_idx], 'b--', linewidth=2, label='3. Planned Outage')
plt.plot(Time_Week2, PU4[h_start_idx:h_end_idx], 'r:', linewidth=2, label='4. Unplanned Outage')

# Highlight the 24h outage window
outage_hours = Time[h_outage_start_idx:h_outage_end_idx]
if len(outage_hours) > 0:
	plt.axvspan(outage_hours[0], outage_hours[-1], color='orange', alpha=0.2, label='24h Thermal Outage Window')

plt.title('Resilience Impact: Unserved Energy During Dunkelflaute Week (Focus on Outage)')
plt.ylabel('Unserved Energy (MW)')
plt.xlabel(f'Hour of the Year (t={h_start_idx+1} to {h_end_idx})')
plt.legend(loc='upper right')
plt.grid(True, linestyle='--')
plt.xlim(Time_Week2[0], Time_Week2[-1])

plt.tight_layout()
plt.show()
'''

# 8.3 Final Metrics Summary Table
Metrics = pd.DataFrame({
	'Scenario': ['Normal (S1)', 'Dunkelflaute (S2)', 'Planned Outage (S3)', 'Unplanned Outage (S4)'],
	'Total Cost (£)': [C1, C2, C3, C4],
	'Total Unserved MWh': [np.sum(PU1), np.sum(PU2), np.sum(PU3), np.sum(PU4)]
})

print('\n----------------------------------------------')
print(' 	 ENERGY SYSTEM RESILIENCE METRICS 	 ')
print('----------------------------------------------')
print(f"Total time horizon: {T} hours (28 days).")
print(f"VoLL (Cost of Shedding): £{C_VoLL:,.0f}/MWh.")
print(Metrics.to_markdown(index=False, floatfmt=".2f"))