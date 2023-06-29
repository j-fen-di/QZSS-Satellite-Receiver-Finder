from geopy.geocoders import Nominatim
from formulas import formulas
import numpy as np

''' NOTE: This is a command line script that takes in the latitude, longitude, and altitude of 4 QZSS
    satellites, as well as measured satellite pseudoranges, and calculates the geographic latitude and longitude
    of the receiver location. This will also tell the physical 'address' of the receiver location using a geolocation
    API
    
    Created by Justin Effendi; June 29, 2023
    
    Simple QZSS Satellite Receiver Finder, v. 0.1'''

# Intro message
print(f"\n\nThe QZSS Satellite Receiver Finder, v. 0.1 (simple version) \nCreated by Justin Effendi, Copyright \xA9 2023 \n\n\
SUMMARY: Input the latitude, longitude, and altitude of each of the four QZSS satellites, as well as \n\
the pseudoranges of the four QZSS satellites. In turn, you will get the geographic latitude and longitude, \n\
as well as the physical 'address', of your receiver location using a geolocation API and some calculations. \n\n")

# ask user for latitude, longitude, and altitude of each satellite
# QZSS-1 (lat, lon, h)
lat_qz1 = input("Enter latitude of 1st QZSS satellite (QZSS-1): ")
lon_qz1 = input("Enter longitude of 1st QZSS satellite (QZSS-1): ")
alt_qz1 = input("Enter altitude of 1st QZSS satellite (QZSS-1): ")

# QZSS-2 (lat, lon, h)
lat_qz2 = input("Enter latitude of 2nd QZSS satellite (QZSS-2): ")
lon_qz2 = input("Enter longitude of 2nd QZSS satellite (QZSS-2): ")
alt_qz2 = input("Enter altitude of 2nd QZSS satellite (QZSS-2): ")

# QZSS-3 (lat, lon, h)
lat_qz3 = input("Enter latitude of 3rd QZSS satellite (QZSS-3): ")
lon_qz3 = input("Enter longitude of 3rd QZSS satellite (QZSS-3): ")
alt_qz3 = input("Enter altitude of 3rd QZSS satellite (QZSS-3): ")

# QZSS-2 (lat, lon, h)
lat_qz4 = input("Enter latitude of 4th QZSS satellite (QZSS-4): ")
lon_qz4 = input("Enter longitude of 4th QZSS satellite (QZSS-4): ")
alt_qz4 = input("Enter altitude of 4th QZSS satellite (QZSS-4): ")

# measured pseudorange of each of the 4 QZSS satellites
psu_mes_1 = input("Enter measured pseudorange of 1st QZSS satellite (QZSS-1): ")
psu_mes_2 = input("Enter measured pseudorange of 2nd QZSS satellite (QZSS-2): ")
psu_mes_3 = input("Enter measured pseudorange of 3rd QZSS satellite (QZSS-3): ")
psu_mes_4 = input("Enter measured pseudorange of 4th QZSS satellite (QZSS-4): ")

# organize input variables
''' initial values for latitude (phi), longitude (lambda), and altitude (h) for the four satellites.
    each subvector represents a satellite (i.e. QZS-1, QZS-2, etc.) '''
QZSS = [[float(lat_qz1), float(lon_qz1), float(alt_qz1)],
        [float(lat_qz2), float(lon_qz2), float(alt_qz2)],
        [float(lat_qz3), float(lon_qz3), float(alt_qz3)],
        [float(lat_qz4), float(lon_qz4), float(alt_qz4)]]

# measured (actual) pseudoranges of 4 QZSS satellites
rhos_act = [float(psu_mes_1), float(psu_mes_2), float(psu_mes_3), float(psu_mes_4)]

# initialize formulas class
func = formulas()

# initialize and allocate space for array of QZSS satellite positions in ECEF coordinates
QZSS_ECEF = np.zeros([4, 3])

# find ECEF position for each satellite and print results
for i in range(0, len(QZSS_ECEF)):
    # calculate position in ECEF frame with to_r_ECEF
    r_ECEF = func.to_r_ECEF(QZSS[i][0], QZSS[i][1], QZSS[i][2])
    # add to QZSS_ECEF array for future usage
    QZSS_ECEF[i][:] = r_ECEF[:]

# calculate expected pseudoranges of the 4 QZSS satellites
rhos_exp = []
for i in range(len(QZSS_ECEF)):
    rhos_exp.append(func.calc_pseudo(QZSS_ECEF[i]))

# initial guess <x_u, y_u, z_u, t_u>
guess_0 = np.zeros([1, 4])

# linearize to user location in ECEF coordinates
receiver_pos = func.trilaterate_to_loc(guess_0, QZSS_ECEF, rhos_exp, rhos_act)

# convert to lat, lon, h
[lat, lon, h] = func.findLatLonH(receiver_pos[0][0], receiver_pos[0][1], receiver_pos[0][2])
receiver_loc = [lat, lon, h]

# print results, gives coordinates of Takesaki Observation Stand in Japan
print("\n\n\nRESULTS:")
print(f"The latitude of the receiver is {receiver_loc[0]:.8g} degrees.")
print(f"The longitude of the receiver is {receiver_loc[1]:.9g} degrees.")
print(f"The altitude of the receiver is {receiver_loc[2]:.9g} km.")