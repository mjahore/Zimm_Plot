#!/opt/conda/bin/python3
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Linear function
def linear (x, m, b):
	return m*x + b

# Output the Zimm plot. We have the slopes and y-intercepts
# stored for each set of data. Let's plot the data first and
# the linear fits, and then we'll plot our extrapolations.
fig, ax  = plt.subplots(dpi=150)
pnt_size = 10

# Instrument parameters
laser_wavelength = 546.0  # nm
refractive_index = 1.502  # Benzene
dn_dc            = 0.108  # Polystyrene (?) in benzene
shift_factor     = 5.0    # For shifting curves on the Zimm plot

# Optical constant
K = 4.0 * (3.14159)**2 * (refractive_index)**2 * (dn_dc)**2 / (laser_wavelength * 10**-7)**4 / (6.02 * 10**23)

print("\nOptical constant K: " + str(K) + " cm^2 mol/g^2\n")

# Number of concentrations per angle.
n_angle = int(sys.argv[2]) 
n_conc  = int(sys.argv[3])

# Input file, should contain columns of concentration, R_theta (angle).
data_in = sys.argv[1]
ls_df   = pd.read_csv(data_in, comment='#')

# We want to convert the concentrations from mg/mL to g/mL.
ls_df["c"] = ls_df["c"] / 1000.0

# The Zimm plot needs Kc/R_theta.
ls_df["Kc_R_theta"] = K*ls_df["c"]/ls_df["R_theta"]

# We want to plot these data as a function of q^2, not angle. So
# convert theta to radians, and compute q^2 (units: nm^-2).
ls_df["angle"] = ls_df["angle"] * 3.14159/180.0
ls_df["angle"] = (4.0 * 3.14159 * refractive_index / laser_wavelength * np.sin(ls_df["angle"]/2.0))**2

# Back up the initial data for plotting purposes.
ls_df_orig = ls_df

# Do the extrapolations. First to c = 0.
c_zero_val = np.zeros(n_angle)
c_zero_slp = np.zeros(n_angle)
for i in range(0, n_angle):
	# Select all common angle data.
	Kc_R_theta = ls_df.loc[i*n_conc:i*n_conc+n_conc-1][["c", "angle", "Kc_R_theta"]]

	# Kc_R_theta now contains Kc_R_theta vs. c for each angle. We fit to extrapolate to
	# c -> 0.
	fit_parms, fit_cov = np.polyfit(Kc_R_theta["c"], Kc_R_theta["Kc_R_theta"], 1, cov="unscaled")

	# Get the current angle as an index.
	angle = int(Kc_R_theta["angle"].iloc[0])

	# Record the y-intercept (Kc/R_theta at c = 0):
	c_zero_slp[i] = fit_parms[0]
	c_zero_val[i] = fit_parms[1]

# Insert the c = 0 values of Kc/R_theta to the main dataframe.
for i in range(0, n_angle):
	start_index = i*n_conc + i
	angle       = ls_df["angle"].iloc[start_index]
	c_zero_row  = pd.DataFrame({"c": [0.00], "angle": [angle], "R_theta": [0.00], "Kc_R_theta": [c_zero_val[i]]}) 
	ls_df       = pd.concat([ls_df.iloc[:start_index], c_zero_row, ls_df.iloc[start_index:]]).reset_index(drop=True)

# There is now one additional concentration!
n_conc = n_conc + 1

# Do the extrapolation to q = 0.
q_zero_val = np.zeros(n_conc)
q_zero_slp = np.zeros(n_conc)
for i in range(0, n_conc):
	q_zero = pd.DataFrame()
	for j in range (0, n_angle):
		# Select all c = 0 data for each angle.
		idx = j*n_conc + i
		Kc_R_theta = ls_df.loc[idx][["c", "angle", "Kc_R_theta"]]

		q_zero_row = pd.DataFrame({"c": [Kc_R_theta["c"]], "angle": [Kc_R_theta["angle"]], "Kc_R_theta": [Kc_R_theta["Kc_R_theta"]]})
		q_zero     = pd.concat([q_zero, q_zero_row])

	# Fit to get y-intercept at q = 0.
	fit_parms, fit_cov = np.polyfit(q_zero["angle"], q_zero["Kc_R_theta"], 1, cov="unscaled")
	
	# Record the value.
	q_zero_slp[i] = fit_parms[0]
	q_zero_val[i] = fit_parms[1]	
	

# Record the zero angle data.
q_zero_data = pd.DataFrame()
for i in range(0, n_conc):
	# Get the concentration this was done at.
	conc        = ls_df["c"].iloc[i]
	q_zero_row  = pd.DataFrame({"c": [conc], "angle": [0.00], "R_theta": [0.00], "Kc_R_theta": [q_zero_val[i]]})
	q_zero_data = pd.concat([q_zero_data, q_zero_row])

# We have another angle now!
n_angle = n_angle + 1

# Add the zero angle data.
ls_df = pd.concat([q_zero_data, ls_df]).reset_index(drop=True)

# Calculate Mw from q = 0, c = 0 data:
q_c_zero = ls_df.iloc[0]
Mw = 1.0 / q_c_zero["Kc_R_theta"]
print("Mw: " + str(math.trunc(Mw)) + " g/mol\n")

# Calculate the radius of gyration from c = 0 slope.
Rg2 = 3.0 * Mw * q_zero_slp[0]
Rg  = np.sqrt(Rg2)
print ("Rg: " + str(math.trunc(Rg*100)/100) + " nm.\n")

# Calculate the z-avg. 2nd virial coefficient.
q_zero_data = pd.DataFrame()
for i in range(0, n_angle):
	q_zero_row  = ls_df.loc[i][["c", "angle", "Kc_R_theta"]]
	q_zero_row  = pd.DataFrame({"c": [q_zero_row["c"]], "angle": [q_zero_row["angle"]], "Kc_R_theta": [q_zero_row["Kc_R_theta"]]})
	q_zero_data = pd.concat([q_zero_data, q_zero_row]).reset_index(drop=True)

# Fit to get slope in the q = 0 condition.
fit_parms, fit_cov = np.polyfit(q_zero_data["c"], q_zero_data["Kc_R_theta"], 1, cov="unscaled")
A2 = fit_parms[0] / 2.0

# We need to do some unit conversion 
print("A2: " + str(math.trunc(1000000*A2)/1000000) + " cm^3 mol/g^2\n")

# Create our scaled x-axis.
ls_df["qc"]      = ls_df["angle"] + shift_factor*ls_df["c"]
ls_df_orig["qc"] = ls_df_orig["angle"] + shift_factor*ls_df_orig["c"]

# This is all data.
ax.scatter(ls_df_orig["qc"], ls_df_orig["Kc_R_theta"], marker='o', label='Data', color='black', s=pnt_size)

# Add the c = 0 extrapolations.
c_zero_data = pd.DataFrame()
for i in range(0, n_angle):
	idx = i*n_conc
	c_zero_row  = ls_df.loc[idx][["c", "angle", "qc", "Kc_R_theta"]]
	c_zero_row  = pd.DataFrame({"c": [c_zero_row["c"]], "angle": [c_zero_row["angle"]], "qc": [c_zero_row["qc"]], "Kc_R_theta": [c_zero_row["Kc_R_theta"]]})
	c_zero_data = pd.concat([c_zero_data, c_zero_row])

print(c_zero_data)
print("\n")

# Plot the extrapolated points.
ax.scatter(c_zero_data["qc"], c_zero_data["Kc_R_theta"], marker='x', label='c = 0 extrapolation', color='red', s=pnt_size)

# Add linear fit.
fit_parms, fit_cov = np.polyfit(c_zero_data["angle"], c_zero_data["Kc_R_theta"], 1, cov="unscaled")
c_zero = ax.plot(c_zero_data["qc"], linear(c_zero_data["angle"], fit_parms[0], fit_parms[1]), color='red', linestyle='dashed')

# Add the q = 0 extrapolations.
q_zero_data = pd.DataFrame()
for i in range(0, n_angle):
	q_zero_row  = ls_df.loc[i][["c", "angle", "qc", "Kc_R_theta"]]
	q_zero_row  = pd.DataFrame({"c": [q_zero_row["c"]], "angle": [q_zero_row["angle"]], "qc": [q_zero_row["qc"]], "Kc_R_theta": [q_zero_row["Kc_R_theta"]]})
	q_zero_data = pd.concat([q_zero_data, q_zero_row])

print(q_zero_data)

# Plot the extrapolated points.
ax.scatter(q_zero_data["qc"], q_zero_data["Kc_R_theta"], marker='^', label='q = 0 extrapolation', color='blue', s=pnt_size)

# Add the linear fits.
fit_parms, fit_cov = np.polyfit(q_zero_data["c"], q_zero_data["Kc_R_theta"], 1, cov="unscaled")
q_zero = ax.plot(q_zero_data["qc"], linear(q_zero_data["c"], fit_parms[0], fit_parms[1]), color='blue', linestyle='dashed')

# Add black lines for the remaining fits.
#
# First the fits wrt concentration:
for i in range(1, n_angle):
	angle_row = ls_df.loc[i*(n_conc):i*(n_conc)+(n_conc-1)]
	fit_parms, fit_cov = np.polyfit(angle_row["c"], angle_row["Kc_R_theta"], 1, cov="unscaled")
	ax.plot(angle_row["qc"], linear(angle_row["c"], fit_parms[0], fit_parms[1]), color='black', linestyle='solid')

# Second, the fits wrt q:
for i in range(1, n_conc):
	c_frame = pd.DataFrame()
	for j in range (1, n_angle):
		idx = j*n_conc + i
		c_row = ls_df.loc[idx]
		c_row = pd.DataFrame({"c": [c_row["c"]], "qc": [c_row["qc"]], "Kc_R_theta": [c_row["Kc_R_theta"]], "angle": [c_row["angle"]]})
		c_frame = pd.concat([c_frame, c_row])

	fit_parms, fit_cov = np.polyfit(c_frame["angle"], c_frame["Kc_R_theta"], 1, cov="unscaled")
	ax.plot(c_frame["qc"], linear(c_frame["angle"], fit_parms[0], fit_parms[1]), color='black', linestyle='solid')

		
# Label axes.
ax.set_ylabel(r'$Kc/R_{\theta}$ (mol/g)')
ax.set_xlabel(r'$q^{2} + ' + str(shift_factor) + r'\times c$ (nm$^{-2}$)')

# Display everything.
plt.legend()
plt.show()

# Print the completed dataframe with c = 0 and q/theta = 0 values.
print("\nComplete dataframe with extrapolated values:\n")
print(ls_df)
print("\n")

# Pat yourself on the back for doing nothing and letting the computer
# do all of the work.
quit

