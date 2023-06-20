import numpy as np

# NOTE: This is a re-implementation of MATLAB code from AE 4361, HW 5, Q2, but in Python :)
# jfendi's code for AE 4361, Homework 5, Question 2 - Spring 2022 (originally MATLAB, now Python)
# Justin Effendi, Dr. Lightsey, AE 4361, February 24, 2022

# PURPOSE: function converts to lat, lon, h to ECEF reference frame
def to_r_ECEF(lat, lon, h):
    r_E = 6371 # radius of Earth, in km
    x = (r_E + h) * np.cos(np.rad2deg(lon)) * np.cos(np.rad2deg(lat))
    y = (r_E + h) * np.sin(np.rad2deg(lon)) * np.cos(np.rad2deg(lat))
    z = (r_E + h) * np.sin(np.rad2deg(lat))
    r_ECEF = [x, y, z]
    return r_ECEF

# PROBLEM 2A
# initial values for latitude (phi), longitude (lambda), and altitude (h)
# for the four satellites
# NOTE: each subvector represents a satellite (i.e. QZS-1, QZS-2, etc.)
QZSS = [[-11.28, 127.27, 34666.42],
        [41.68, 137.46, 38886.86],
        [0.05, 127.02, 35790.31],
        [-11.12, 149.14, 34588.99]]

# initialize and allocate space for array of QZSS satellite positions
QZSS_ECEF = np.zeros([4, 3])

# print '2a)'
print("2a)")

# find ECEF position for each satellite and print results (2a)
for i in range(0, len(QZSS_ECEF)):
    # calculate position in ECEF frame with to_r_ECEF
    r_ECEF = to_r_ECEF(QZSS[i][0], QZSS[i][1], QZSS[i][2])
    # add to QZSS_ECEF array for future usage
    QZSS_ECEF[i][:] = r_ECEF[:]
    # print results
    res2a = f"r of the QZS-{i} satellite is <{r_ECEF[0]:.3f}, {r_ECEF[1]:.3f}, {r_ECEF[2]:.3f}> km.\n"
    print(res2a)

# PROBLEM 2B
# print '2b)'
print('2b)')
r_origin_ECEF = [0, 0, 0]
c = 299792.000; # in km/s
t_u = 0; # in seconds (1 day)
for i in range(0, len(QZSS_ECEF)):
    rho_hat = np.sqrt((QZSS_ECEF[i][0])**2 + (QZSS_ECEF[i][1])**2 + (QZSS_ECEF[i][2])**2) + c*t_u
    # print results
    res2b = f"pseudorange of the QZS-{i} satellite is {rho_hat:.3f} km.\n"
    print(res2b)