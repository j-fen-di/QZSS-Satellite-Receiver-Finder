import numpy as np

'''NOTE: These are the formulas, algorithms, and convergence methods that assist in calculating
the geographic latitude and longitude of a receiver's location based on the latitude, longitude,
altitude, and pseudorange of each of the 4 QZSS satellites.'''

# CONSTANTS
c = 299792.000 # speed of light (in km/s)
r_E = 6371 # radius of Earth (in km)
t_u = 0 # expected clock receiver offset (in seconds - assume there is none)

class formulas:
    # def __init__(self, QZSS_locs, QZSS_psu):
    #     self.QZSS_locs = QZSS_locs
    #     self.QZSS_psu = QZSS_psu

    # PURPOSE: function converts to lat, lon, h to ECEF reference frame
    def to_r_ECEF(self, lat, lon, h):
        x = (r_E + h) * np.cos(np.deg2rad(lon)) * np.cos(np.deg2rad(lat))
        y = (r_E + h) * np.sin(np.deg2rad(lon)) * np.cos(np.deg2rad(lat))
        z = (r_E + h) * np.sin(np.deg2rad(lat))
        r_ECEF = [x, y, z]
        return r_ECEF
    
    # PURPOSE: calculate pseudorange of each QZSS satellite based on ECEF lat, lon, and altitude
    def calc_pseudo(self, sat):
        rho_hat = np.sqrt((sat[0])**2 + (sat[1])**2 + (sat[2])**2) + c*t_u
        return rho_hat
    
    def trilaterate_to_loc(self, r_user, r_sats, pseudo_exp, pseudo_act):
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
        if formulas.isConvergence(self, pseudo_guess, pseudo_act) == 1:
            r_user_f = r_user
        else:
            r_user_f = formulas.trilaterate_to_loc(self, r_user, r_sats, pseudo_guess, pseudo_act) # recurse until condition met

        return r_user_f

    def isConvergence(self, a, b):
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

    def findLatLonH(self, x, y, z):
        # converts ECEF vector to lat, lon, h on Earth surface
        # needed constants
        r_E = 6371 # radius of Earth, in km
        
        # solving for lon, lat, and h using algebra
        lon = np.rad2deg(np.arctan2(y, x)) # calculating latitude (in degrees)
        lat = np.rad2deg(np.arctan2(z * np.sin(np.deg2rad(lon)), y)) # calculating longitude (in degrees)
        h = z/np.sin(np.deg2rad(lat)) - r_E # calculating altitude (in meters)

        return lat, lon, h