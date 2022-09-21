from math import *
from numpy import *
import scipy as scp

#
# --------------------------------------------------------------------- 
# |                                                                     |
# |	geodetic.cc -  a collection of geodetic functions                   |
# |	Jim Leven  - Dec 99                                                 |
# |                                                                     |
# | originally from:                                                    |
# | http://wegener.mechanik.tu-darmstadt.de/GMT-Help/Archiv/att-8710/Geodetic_py |                                                                   |
# |                                                                     |
# --------------------------------------------------------------------- 
# 
# 
# ----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual	    |
# | 								                                    |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html	        		|
# | 								                                    |
# | This page last updated 11 May 1999 	                				|
# | 								                                    |
# | Computations on the Ellipsoid	                    				|
# | 								                                    |
# | There are a number of formulae that are available           		|
# | to calculate accurate geodetic positions, 		            		|
# | azimuths and distances on the ellipsoid.			                |
# | 								                                    |
# | Vincenty's formulae (Vincenty, 1975) may be used 		            |
# | for lines ranging from a few cm to nearly 20,000 km, 	            |
# | with millimetre accuracy. 					                        |
# | The formulae have been extensively tested 		                    |
# | for the Australian region, by comparison with results       		|
# | from other formulae (Rainsford, 1955 & Sodano, 1965). 	            |
# |								                                        |
# | * Inverse problem: azimuth and distance from known 	        		|
# |			latitudes and longitudes 			                        |
# | * Direct problem: Latitude and longitude from known 	            |
# |			position, azimuth and distance. 		                    |
# | * Sample data 						                                |
# | * Excel spreadsheet 			                            		|
# | 								                                    |
# | Vincenty's Inverse formulae				                    		|
# | Given: latitude and longitude of two points                 		|
# |			(phi1, lembda1 and phi2, lembda2), 	                        |
# | Calculate: the ellipsoidal distance (s) and 	            		|
# | forward and reverse azimuths between the points (alpha12, alpha21).	|
# |									                                    |
# ---------------------------------------------------------------------- 

def vinc_dist(  phi1,  lembda1,  phi2,  lembda2 ) :
        """ 
        Returns the distance between two geographic points on the ellipsoid
        and the forward and reverse azimuths between these points.
        lats, longs and azimuths are in decimal degrees, distance in metres 

        Returns ( s, alpha12,  alpha21 ) as a tuple
        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        if (abs( phi2 - phi1 ) < 1e-8) and ( abs( lembda2 - lembda1) < 1e-8 ) :
                return 0.0, 0.0, 0.0

        piD4   = math.atan( 1.0 )
        two_pi = piD4 * 8.0

        phi1    = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0		# unfortunately lambda is a key word!
        phi2    = phi2 * piD4 / 45.0
        lembda2 = lembda2 * piD4 / 45.0

        b = a * (1.0 - f)

        TanU1 = (1-f) * tan( phi1 )
        TanU2 = (1-f) * tan( phi2 )

        U1 = atan(TanU1)
        U2 = atan(TanU2)

        lembda = lembda2 - lembda1
        last_lembda = -4000000.0		# an impossibe value
        omega = lembda

        # Iterate the following equations, 
        #  until there is no significant change in lembda 

        while ( last_lembda < -3000000.0 or lembda != 0 and abs( (last_lembda - lembda)/lembda) > 1.0e-9 ) :

                sqr_sin_sigma = pow( cos(U2) * sin(lembda), 2) + \
                        pow( (cos(U1) * sin(U2) - \
                        sin(U1) *  cos(U2) * cos(lembda) ), 2 )

                Sin_sigma = sqrt( sqr_sin_sigma )

                Cos_sigma = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(lembda)
        
                sigma = atan2( Sin_sigma, Cos_sigma )

                Sin_alpha = cos(U1) * cos(U2) * sin(lembda) / sin(sigma)
                alpha = asin( Sin_alpha )

                Cos2sigma_m = cos(sigma) - (2 * sin(U1) * sin(U2) / pow(cos(alpha), 2) )

                C = (f/16) * pow(cos(alpha), 2) * (4 + f * (4 - 3 * pow(cos(alpha), 2)))

                last_lembda = lembda

                lembda = omega + (1-C) * f * sin(alpha) * (sigma + C * sin(sigma) * \
                        (Cos2sigma_m + C * cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

        u2 = pow(cos(alpha),2) * (a*a-b*b) / (b*b)

        A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

        B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))

        delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * \
                (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) ) - \
                (B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * \
                (-3 + 4 * pow(Cos2sigma_m,2 ) )))

        s = b * A * (sigma - delta_sigma)

        alpha12 = atan2( (cos(U2) * sin(lembda)), \
                (cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lembda)))

        alpha21 = atan2( (cos(U1) * sin(lembda)), \
                (-sin(U1) * cos(U2) + cos(U1) * sin(U2) * cos(lembda)))

        if ( alpha12 < 0.0 ) : 
                alpha12 =  alpha12 + two_pi
        if ( alpha12 > two_pi ) : 
                alpha12 = alpha12 - two_pi

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) : 
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) : 
                alpha21 = alpha21 - two_pi

        alpha12    = alpha12    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4
        return s, alpha12,  alpha21 

   # END of Vincenty's Inverse formulae 


#-------------------------------------------------------------------------------
# Vincenty's Direct formulae							|
# Given: latitude and longitude of a point (phi1, lembda1) and 			|
# the geodetic azimuth (alpha12) 						|
# and ellipsoidal distance in metres (s) to a second point,			|
# 										|
# Calculate: the latitude and longitude of the second point (phi2, lembda2) 	|
# and the reverse azimuth (alpha21).						|
# 										|
#-------------------------------------------------------------------------------

def  vinc_pt( phi1, lembda1, alpha12, s ) :
        """

        Returns the lat and long of projected point and reverse azimuth
        given a reference point and a distance and azimuth to project.
        lats, longs and azimuths are passed in decimal degrees

        Returns ( phi2,  lambda2,  alpha21 ) as a tuple 

        """
 
        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        piD4 = atan( 1.0 )
        two_pi = piD4 * 8.0

        phi1    = phi1    * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0
        alpha12 = alpha12 * piD4 / 45.0
        if ( alpha12 < 0.0 ) : 
                alpha12 = alpha12 + two_pi
        if ( alpha12 > two_pi ) : 
                alpha12 = alpha12 - two_pi

        b = a * (1.0 - f)

        TanU1 = (1-f) * tan(phi1)
        U1 = atan( TanU1 )
        sigma1 = atan2( TanU1, cos(alpha12) )
        Sinalpha = cos(U1) * sin(alpha12)
        cosalpha_sq = 1.0 - Sinalpha * Sinalpha

        u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
        A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
                (320 - 175 * u2) ) )
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

        # Starting with the approximation
        sigma = (s / (b * A))

        last_sigma = 2.0 * sigma + 2.0	# something impossible

        # Iterate the following three equations 
        #  until there is no significant change in sigma 

        # two_sigma_m , delta_sigma
        while ( abs( (last_sigma - sigma) / sigma) > 1.0e-9 ) :
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = B * sin(sigma) * ( cos(two_sigma_m) \
                        + (B/4) * (cos(sigma) * \
                        (-1 + 2 * pow( cos(two_sigma_m), 2 ) -  \
                        (B/6) * cos(two_sigma_m) * \
                        (-3 + 4 * pow(sin(sigma), 2 )) *  \
                        (-3 + 4 * pow( cos (two_sigma_m), 2 ))))) \

                last_sigma = sigma
                sigma = (s / (b * A)) + delta_sigma

        phi2 = atan2 ( (sin(U1) * cos(sigma) + cos(U1) * sin(sigma) * cos(alpha12) ), \
                ((1-f) * sqrt( pow(Sinalpha, 2) +  \
                pow(sin(U1) * sin(sigma) - cos(U1) * cos(sigma) * cos(alpha12), 2))))

        lembda = atan2( (sin(sigma) * sin(alpha12 )), (cos(U1) * cos(sigma) -  \
                sin(U1) *  sin(sigma) * cos(alpha12)))

        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

        omega = lembda - (1-C) * f * Sinalpha *  \
                (sigma + C * sin(sigma) * (cos(two_sigma_m) + \
                C * cos(sigma) * (-1 + 2 * pow(cos(two_sigma_m),2) )))

        lembda2 = lembda1 + omega

        alpha21 = atan2 ( Sinalpha, (-sin(U1) * sin(sigma) +  \
                cos(U1) * cos(sigma) * cos(alpha12)))

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) :
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) :
                alpha21 = alpha21 - two_pi

        phi2       = phi2       * 45.0 / piD4
        lembda2    = lembda2    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4

        return phi2,  lembda2,  alpha21 

  # END of Vincenty's Direct formulae

#--------------------------------------------------------------------------
# Notes: 
# 
# * "The inverse formulae may give no solution over a line 
# 	between two nearly antipodal points. This will occur when 
# 	lembda ... is greater than pi in absolute value". (Vincenty, 1975)
#  
# * In Vincenty (1975) L is used for the difference in longitude, 
# 	however for consistency with other formulae in this Manual, 
# 	omega is used here. 
# 
# * Variables specific to Vincenty's formulae are shown below, 
# 	others common throughout the manual are shown in the Glossary. 
# 
# 
# alpha = Azimuth of the geodesic at the equator
# U = Reduced latitude
# lembda = Difference in longitude on an auxiliary sphere (lembda1 & lembda2 
# 		are the geodetic longitudes of points 1 & 2)
# sigma = Angular distance on a sphere, from point 1 to point 2
# sigma1 = Angular distance on a sphere, from the equator to point 1
# sigma2 = Angular distance on a sphere, from the equator to point 2
# sigma_m = Angular distance on a sphere, from the equator to the 
# 		midpoint of the line from point 1 to point 2
# u, A, B, C = Internal variables
# 
# 

#*******************************************************************
