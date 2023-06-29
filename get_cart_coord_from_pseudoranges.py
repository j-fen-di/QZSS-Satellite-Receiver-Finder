import numpy as np

# NOTE: This is a re-implementation of MATLAB code from AE 4361, HW 5, Q3, but in Python :)
# jfendi's MATLAB code for AE 4361, Homework 5, Question 3 - Spring 2022
# Justin Effendi, Dr. Lightsey, AE 4361, February 24, 2022

# format long (Do something about this lol)

# FUNCTIONS
def trilaterate_to_loc(r_user, r_sats, pseudo_exp, pseudo_act):
    # initial values needed
    c = 299792.000 # speed of light, in km/s
    
    # measured - expected pseudoranges
    delta_rho = np.subtract(pseudo_act, pseudo_exp)
    
    # initialize matrix H
    H = [[0, 0, 0, c],
        [0, 0, 0, c],
        [0, 0, 0, c],
        [0, 0, 0, c]]
    
    # edit H to proper form
    for i in range(4):
        r_i_hat = np.sqrt((r_sats[i][0] - r_user[0][0])**2 + (r_sats[i][1] - r_user[0][1])**2 + (r_sats[i][2] - r_user[0][2])**2)
        for j in range(3):
            a = -(r_sats[i][j] - r_user[0][j]) / r_i_hat
            H[i][j] = a
    
    # calculate correction term (inverse H times transposed delta_rho)
    delta_x_u = np.matmul(np.linalg.inv(H), np.transpose(delta_rho))
    delta_x_u = np.transpose(delta_x_u); # transpose delta_x_u
    
    # update r_user with correction term
    r_user = np.add(r_user, delta_x_u)
    
    # come up with pseudo_guess
    pseudo_guess = np.zeros([1, 4]) # initialize pseudo_guess
    for i in range(4):
        pseudo_guess[0][i] = np.sqrt((r_sats[i][0] - r_user[0][0])**2 + (r_sats[i][1] - r_user[0][1])**2 + 
                                  (r_sats[i][2] - r_user[0][2])**2) + c*r_user[0][3]
    
    # set final answer if converges, recurse otherwise
    if isConvergence(pseudo_guess, pseudo_act) == 1:
        r_user_f = r_user
    else:
        r_user_f = trilaterate_to_loc(r_user, r_sats, pseudo_guess, pseudo_act) # recurse until condition met

    return r_user_f

def isConvergence(a, b):
    # function checking if convergence criteria is met
    x = 0 # pseudo boolean variable
    switchy = 0 # true if not converge, false if converges
    
    # compare each element of the two vectors
    for i in range(len(a)):
        if (a[i] - b[i]).all() > 1e-4:
            switchy = 1
    
    # return true if convergence
    if switchy == 0:
        x = 1

    return x

def findLatLonH(x, y, z):
    # converts ECEF vector to lat, lon, h on Earth surface
    # needed constants
    r_E = 6371 # radius of Earth, in km
    
    # solving for lon, lat, and h using algebra
    lon = np.rad2deg(np.arctan2(y, x)) # calculating latitude (in degrees)
    lat = np.rad2deg(np.arctan2(z * np.sin(np.deg2rad(lon)), y)) # calculating longitude (in degrees)
    h = z/np.sin(np.deg2rad(lat)) - r_E # calculating altitude (in meters)

    return lat, lon, h

# measured pseudoranges for QZS 1, 2, 3, and 4 (in km)
rhos_act = [36536.1926, 39061.5413, 36817.8847, 36735.0892]

# expected pseudoranges from problem 2b (in km)
rhos_exp = [41037.420, 45257.860, 42161.310, 40959.990]

# satellite conditions for satellites 1, 2, 3, 4 respectively of QZSS
r_QZSS = [[-24371.056, 32026.357, -8027.077],
    [-24905.318, 22853.526, 30095.105],
    [-25385.052, 33662.647, 36.793],
    [-34500.873, 20615.646, -7899.728]]

# initial guess <x_u, y_u, z_u, t_u>
guess_0 = np.zeros([1, 4])

# linearize to user location
receiver_pos = trilaterate_to_loc(guess_0, r_QZSS, rhos_exp, rhos_act)
# print results
print("3)")
res3_1 = f"The location of the receiver in ECEF coordinates are <{receiver_pos[0][0]:.4f}, {receiver_pos[0][1]:.4f}, {receiver_pos[0][2]:.4f}> km."
res3_2 = f"The receiver clock offset is {receiver_pos[0][3]:.10f} seconds.\n\n"
print(res3_1)
print(res3_2)

# BONUS PROBLEM
# convert to lat, lon, h
[lat, lon, h] = findLatLonH(receiver_pos[0][0], receiver_pos[0][1], receiver_pos[0][2])
receiver_loc = [lat, lon, h]
# print results, gives coordinates of Takesaki Observation Stand in Japan
print("Bonus problem:")
print(f"The latitude of the receiver is {receiver_loc[0]:.8g} degrees.")
print(f"The longitude of the receiver is {receiver_loc[1]:.9g} degrees.")
print(f"The altitude of the receiver is {receiver_loc[2]:.9g} km.")