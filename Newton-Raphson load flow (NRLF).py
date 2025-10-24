import pandas as pd
import numpy as np
import math
from New_Ybus import Y  # Your Y-bus function

# --- Load data ---

data1 = pd.read_excel("E:\\folder_content_of_drive_E\\Downloads_DriveC_folders\\Question banks and notes\\M.tech_IITKGP\\1st_SEM\\Power_&_Energy_Lab\\Python_LAB\\Assignment2\\Assignment_2.xlsx", sheet_name="Bus_Data")
line_data = pd.read_excel("E:/folder_content_of_drive_E/Downloads_DriveC_folders/Question banks and notes/M.tech_IITKGP/1st_SEM/Power_&_Energy_Lab/Python_LAB/Assignment2/Assignment_2.xlsx", sheet_name="Line_Data")

# Y-bus matrix
ybus = Y(data1, line_data)
xr = np.array(ybus, dtype=complex)

# --- Extract bus data ---
busdatanr = data1.to_numpy()
n = len(busdatanr)

Bus_No = np.array(data1["Bus_No."], dtype=int)
Bus_type = np.array(data1["Bus_Type"], dtype=int)
Bus_Voltage = np.array(data1["Bus_Voltage_Mag(pu)"], dtype=float)
V_angle = np.array(data1["Bus_Voltage_Ang(pu)"], dtype=float)
V_angle = np.deg2rad(V_angle)  # convert to radians
P_sp = (data1["P_g(MW)"] - data1["P_L(MW)"]).to_numpy() / 100
Q_sp = (data1["Q_g(MVAr)"] - data1["Q_L(MVAr)"]).to_numpy() / 100


# --- Identify bus types ---
slack_list = Bus_No[Bus_type == 0].tolist()
pv_list = Bus_No[Bus_type == 1].tolist()
pq_list = Bus_No[Bus_type == 2].tolist()

p_mismatch_buses = sorted(pv_list + pq_list)
q_mismatch_buses = sorted(pq_list)

# --- Iteration parameters ---
tolerance = 1e-6
max_iter = 1
iteration = 0
DeltaP = np.zeros(n)
DeltaQ = np.zeros(n)


# --- Main Newton-Raphson Loop ---
while iteration < max_iter:
    P_calc = np.zeros(n)
    Q_calc = np.zeros(n)
    for i in range(n):
        for j in range(n):
            Y_ij = xr[i, j]
            P_calc[i] += Bus_Voltage[i]*Bus_Voltage[j]*abs(Y_ij)*np.cos(np.angle(Y_ij) + V_angle[j] - V_angle[i])
            Q_calc[i] += -Bus_Voltage[i]*Bus_Voltage[j]*abs(Y_ij)*np.sin(np.angle(Y_ij) + V_angle[j] - V_angle[i])
    
    # 2. Compute mismatches
    delta_P = np.array([P_sp[b-1] - P_calc[b-1] for b in p_mismatch_buses])
    delta_Q = np.array([Q_sp[b-1] - Q_calc[b-1] for b in q_mismatch_buses])
    A = np.hstack((delta_P, delta_Q))
    
    
    max_mismatch = np.max(np.abs(A))

    print(f"Iteration {iteration+1}: Max mismatch = {max_mismatch:.6f}")
    
    if max_mismatch < tolerance:
        print(f"Converged in {iteration + 1} iterations!")

        #print(f"\n✅ Converged in {iteration} iterations!")
        break

    # 3. Form Jacobian
    # Extract indices
    NonSlack_indices = np.where(Bus_type != 0)[0]
    PQ_indices = np.where(Bus_type == 2)[0]
    n_non = len(NonSlack_indices)
    n_pq = len(PQ_indices)
    
    G = xr.real
    B = xr.imag
    
    # Initialize Jacobian submatrices
    J11 = np.zeros((n_non, n_non))
    J12 = np.zeros((n_non, n_pq))
    J21 = np.zeros((n_pq, n_non))
    J22 = np.zeros((n_pq, n_pq))
    
    # --- J11 ---
    for ii, i in enumerate(NonSlack_indices):
        for jj, j in enumerate(NonSlack_indices):
            if i == j:
                J11[ii, jj] = 0
                for k in range(n):
                    if k != i:
                        J11[ii, jj] += Bus_Voltage[i]*Bus_Voltage[k]*abs(xr[i,k])*np.sin(np.angle(xr[i,k]) + V_angle[k] - V_angle[i])
            else:
                J11[ii, jj] = -Bus_Voltage[i]*Bus_Voltage[j]*abs(xr[i,j])*np.sin(np.angle(xr[i,j]) + V_angle[j] - V_angle[i])
    
    # --- J12 ---
    for ii, i in enumerate(NonSlack_indices):
        for jj, j in enumerate(PQ_indices):
            if i == j:
                sum_term = 0
                for k in range(n):
                    if k != i:
                        sum_term += Bus_Voltage[k]*abs(xr[i,k])*np.cos(np.angle(xr[i,k]) + V_angle[k] - V_angle[i])
                J12[ii, jj] = 2*Bus_Voltage[i]*G[i,i] + sum_term
            else:
                J12[ii, jj] = Bus_Voltage[i]*abs(xr[i,j])*np.cos(np.angle(xr[i,j]) + V_angle[j] - V_angle[i])
    
    # --- J21 ---
    for ii, i in enumerate(PQ_indices):
        for jj, j in enumerate(NonSlack_indices):
            if i == j:
                sum_term = 0
                for k in range(n):
                    if k != i:
                        sum_term += Bus_Voltage[i]*Bus_Voltage[k]*abs(xr[i,k])*np.cos(np.angle(xr[i,k]) + V_angle[k] - V_angle[i])
                J21[ii, jj] = sum_term
            else:
                J21[ii, jj] = -Bus_Voltage[i]*Bus_Voltage[j]*abs(xr[i,j])*np.cos(np.angle(xr[i,j]) + V_angle[j] - V_angle[i])
    
    # --- J22 ---
    for ii, i in enumerate(PQ_indices):
        for jj, j in enumerate(PQ_indices):
            if i == j:
                sum_term = 0
                for k in range(n):
                    if k != i:
                        sum_term += Bus_Voltage[k]*abs(xr[i,k])*np.sin(np.angle(xr[i,k]) + V_angle[k] - V_angle[i])
                J22[ii, jj] = -2*Bus_Voltage[i]*B[i,i] - sum_term
            else:
                J22[ii, jj] = -Bus_Voltage[i]*abs(xr[i,j])*np.sin(np.angle(xr[i,j]) + V_angle[j] - V_angle[i])
    
    # Combine Jacobian
    Jacobian = np.block([[J11, J12],
                         [J21, J22]])
    
    # 4. Solve for correction factors
    correction = np.linalg.solve(Jacobian, A)
    print(correction)
    
    # 5. Update voltage angles and magnitudes
    for idx, bus in enumerate(p_mismatch_buses):
        V_angle[bus-1] += correction[idx]
    for idx, bus in enumerate(q_mismatch_buses):
        Bus_Voltage[bus-1] += correction[len(p_mismatch_buses) + idx]
    
    iteration += 1

# --- Display final voltage magnitudes and angles ---
print("\nFinal Bus Voltages (pu):", Bus_Voltage)
print("Final Bus Angles (deg):", np.rad2deg(V_angle))



def calc_power(V, delta, Ybus):
    n = len(V)
    G = xr.real
    B = xr.imag

    P = np.zeros(n)
    Q = np.zeros(n)

    for i in range(n):
        for j in range(n):
            P[i] += V[i] * V[j] * (G[i, j] * np.cos(delta[i] - delta[j]) +
                                   B[i, j] * np.sin(delta[i] - delta[j]))
            Q[i] += V[i] * V[j] * (G[i, j] * np.sin(delta[i] - delta[j]) -
                                   B[i, j] * np.cos(delta[i] - delta[j]))
    return P, Q
P_final, Q_final = calc_power(Bus_Voltage, V_angle, xr)

print("\nFinal Bus Power Injections:")
for i in range(n):
    print(f"Bus {Bus_No[i]}:  P = {P_final[i]:.6f} pu,   Q = {Q_final[i]:.6f} pu")




