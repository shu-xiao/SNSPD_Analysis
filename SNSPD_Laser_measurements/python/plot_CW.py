#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

Voltage = [0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
pcr_mean = [60.000, 90.000, 1135.000, 1955.000, 2475.000, 4045.000, 41445.000, 100, 50]
pcr_error = [17.720, 30.315, 151.304, 179.038, 212.058, 316.077, 2549.116, 40, 15]

plt.errorbar(Voltage, pcr_mean, yerr=pcr_error, fmt='o', capsize=5)
plt.xlabel('Bias Voltage (V)', fontsize=15)
plt.ylabel('PCR Mean [Hz]', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(True)
plt.grid(True, which='minor', linestyle=':', linewidth='0.5')  # Minor grid
plt.minorticks_on()  # Enable minor ticks
info=r'$T=4.6K,\quad I_{Laser}=100uW,\quad 532nm\, CW\, laser$'
plt.text(1.9, 70000, info, fontsize=13, horizontalalignment='right', verticalalignment='bottom')
plt.tight_layout()
plt.savefig(f"Counts.png")
plt.show()

#Data Set 1
x1 = [0.5, 1.0, 2.0, 4.0, 6.0, 10.0, 16.0, 32.0, 128.0, 150.0]
y1 = [0.00465, 0.00606, 0.0089, 0.01595, 0.01718, 0.02409, 0.0306, 0.04077, 0.04449, 0.04918]
y1_error = [0.00242, 0.00142, 0.00137, 0.00182, 0.0029, 0.00319, 0.00339, 0.00263, 0.00305, 0.00315]

# Data Set 2
x2 = [1.0, 4.0, 8.0, 16.0, 32.0, 64.0]
y2 = [0.05717, 0.05024, 0.05208, 0.0537, 0.0655, 0.07316]
y2_error = [0.00806, 0.01676, 0.01384, 0.01667, 0.01776, 0.00909]

plt.errorbar(x1, y1, yerr=y1_error, fmt='o', label='2V')
plt.errorbar(x2, y2, yerr=y2_error, fmt='o', label='1.7V')

plt.xlabel(r'Laser Intensity [$\mu$W]', fontsize=15)
plt.ylabel(r'Signal Amplitude Mean ($\overline{A}_{signal}$) [V]', fontsize=15)
legend = plt.legend(fontsize='large', loc="lower right")
legend.set_title('Bias Voltage', prop={'size': 'large'})
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(True)
info=r'$T=4.6K,\quad 532nm\, pulsed\, laser$'
plt.text(150, 0.088, info, fontsize=13, horizontalalignment='right', verticalalignment='bottom')
plt.tight_layout()
plt.savefig(f"Compare_Bias_Voltage.png")
plt.show()


####################
####################
####################
####################

x = [0.5, 1, 2]
y = [4.68, 18.875, 56.12]

# Fit a quadratic curve (2nd degree polynomial) to the data
coefficients = np.polyfit(x, y, 2)
curve = np.poly1d(coefficients)

# Generate x values for the curve
x_curve = np.linspace(min(x), max(x), 100)

# Evaluate the curve at the x values
y_curve = curve(x_curve)

# Plot the data points and the quadratic curve
plt.plot(x, y, 'o', label='Data')
plt.plot(x_curve, y_curve, label='Quadratic Fit')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Quadratic Fit of Data Points')

# Add fit equation as text annotation
fit_equation = f'Fit: y = {coefficients[0]:.2f}x^2 + {coefficients[1]:.2f}x + {coefficients[2]:.2f}'
plt.annotate(fit_equation, xy=(0.05, 0.9), xycoords='axes fraction')

plt.legend(loc="lower right")
plt.grid(True)
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Data
x = np.array([0.0, 0.5, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.82, 1.84, 1.86, 1.88, 1.9, 1.95, 2.0])
y = np.array([0.09, 0.075, 0.085, 0.098, 0.116, 0.127, 0.164, 0.23, 0.279, 0.433, 0.651, 0.716, 0.937, 0.976, 0.982, 0.975, 0.994, 0.993])

# Fit with linear spline
f = interp1d(x, y, kind='linear')

# Generate finer x-values for evaluating the fit
x_fine = np.linspace(x.min(), x.max(), 1000)

# Evaluate the fit at finer resolution
y_fine = f(x_fine)

# Find the first index where the fit reaches 0.9
idx = np.where(y_fine >= 0.9)[0][0]
x_90 = x_fine[idx]
y_90 = y_fine[idx]

# Plot the data and the fit
plt.plot(x, y, 'bo', label='Data')
plt.plot(x_fine, y_fine, 'r-', label='Fit')

# Plot the point where fit reaches 0.9
plt.plot(x_90, y_90, 'go', label='0.9 Point')
plt.axvline(x=x_90, color='g', linestyle='--')

# Set labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Linear Spline Fit')

# Add legend
plt.legend()

# Show the plot
plt.show()
