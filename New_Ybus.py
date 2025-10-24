import pandas as pd
import numpy as np

# Read Excel sheets
bus_data = pd.read_excel("E:/folder_content_of_drive_E/Downloads_DriveC_folders/Question banks and notes/M.tech_IITKGP/1st_SEM/Power_&_Energy_Lab/Python_LAB/Assignment2/Assignment_2.xlsx", sheet_name="Bus_Data")
line_data = pd.read_excel("E:/folder_content_of_drive_E/Downloads_DriveC_folders/Question banks and notes/M.tech_IITKGP/1st_SEM/Power_&_Energy_Lab/Python_LAB/Assignment2/Assignment_2.xlsx", sheet_name="Line_Data")


def Y(bus_data, line_data):
    # Clean column names
    bus_data.columns = bus_data.columns.str.strip()
    line_data.columns = line_data.columns.str.strip()

    # Number of buses = highest bus number found in line data
    num_buses = max(line_data['From Bus'].max(), line_data['To Bus'].max())

    # Initialize Ybus matrix
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    Y_L = np.zeros((num_buses, num_buses), dtype=complex)
    
    # Build Ybus
    for _, row in line_data.iterrows():
        from_bus  = int(row['From Bus'])
        to_bus    = int(row['To Bus'])
        r         = row['R']
        x         = row['X']
        ONR=row['Tap Setting']
        ONR = ONR.astype(float)
        # Bus_No = len(from_bus)
        # P_L = row('P_L(MW)')
        # Q_L = row('Q_L(MVAr)')
        # Bus_Voltage = row('Bus_Voltage_Mag(pu)')
        
        
        b_total   = row['Half Line Charging susceptance (p.u.)']
        

        p = from_bus - 1   # index for From bus
        q = to_bus - 1     # index for To bus

        z_series = r + 1j * x
        y_series = 1 / z_series
        y_shunt  = 1j * (b_total / 2.0)
        
        # for x in range(len(Bus_No)):
        #     if (P_L[x]+1j*Q_L[x]) != 0:
        #         Y_L[x,x] = 1/(100 * ((Bus_Voltage[x])**2)/(P_L[x]-1j*Q_L[x]))
                
        # Ybus += Y_L
                

        

        if ONR == 1.0:
            Ybus[p, q] -= y_series
            Ybus[q, p] -= y_series
            Ybus[p, p] += y_series + y_shunt
            Ybus[q, q] += y_series + y_shunt
        else:
            a = ONR
            Ybus[p, p] += y_series / (a**2) + y_shunt
            Ybus[q, q] += y_series + y_shunt
            Ybus[p, q] -= y_series / a
            Ybus[q, p] -= y_series / a

    return Ybus


    


# -------------------------
# # MAIN EXECUTION
# # -------------------------
# if __name__ == "__main__":
#     Ybus = Ybus_c(bus_data, line_data)

#     # Print raw complex matrix
#     print("\nYbus Matrix (complex form):\n", Ybus)

#     # Pretty printing (real + j*imag form like textbooks)
#     print("\nYbus Matrix (separated real & imaginary parts):\n")
#     for i in range(Ybus.shape[0]):
#         row_str = ""
#         for j in range(Ybus.shape[1]):
#             val = Ybus[i, j]
#             row_str += f"{val.real:.4f} + j{val.imag:.4f}\t"
#         print(row_str)
        
      
