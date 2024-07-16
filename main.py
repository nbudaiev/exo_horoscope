#from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import astropy.units as u
from geopy import geocoders  
from geopy.geocoders import Nominatim
from datetime import datetime

def birth_lat_log(citystate):
    geolocator = Nominatim(user_agent='moeur')
    location = geolocator.geocode(citystate)
    birth_lat, birth_log = location[1][0], location[1][1]
    return birth_lat, birth_log

class User(object):
    """
    User class
    """

    def __init__(self, latitude, longitude, time):
        """
        
        """
        if latitude.unit != u.deg:
            raise ValueError('Latitude must be in degrees')
        if longitude.unit != u.deg:
            raise ValueError('Longitude must be in degrees')

        if latitude.value < -90 or latitude.value > 90:
            raise ValueError('Latitude must be between -90 and 90')
        if longitude.value < -180 or longitude.value > 180:
            raise ValueError('Longitude must be between -180 and 180')
        
        self.latitude = latitude
        self.longitude = longitude
        
        self.time = time
        # implement the time using year, month, day,... version (similar to below.) Change the __init__ parameters
        #self.time = Time(datetime(year, month, day, hours, minute, second))
        

        self.ra, self.dec = self.get_zenith()




    def get_zenith(self):
        """
        Get zenith coordinates
        """

        location = EarthLocation(lat=self.latitude, lon=self.longitude)
        zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=AltAz(obstime=self.time, location=location))
        zenith_radec = zenith.transform_to('icrs')



        return zenith_radec.ra, zenith_radec.dec


#test_object = User(48.85*u.deg, 2.35*u.deg, Time.now())
#print(test_object.latitude)
