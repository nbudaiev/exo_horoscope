#from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import astropy.units as u



class User(object):
    """
    User class
    """

    def __init__(self, latitude, longitude, time):
        """
        
        """
        # if latitude.unit != 'u.deg':
        #     raise ValueError('Latitude must be in degrees')
        # if longitude.unit != 'u.deg':
        #     raise ValueError('Longitude must be in degrees')

        # if latitude.value < -90 or latitude.value > 90:
        #     raise ValueError('Latitude must be between -90 and 90')
        # if longitude.value < -180 or longitude.value > 180:
        #     raise ValueError('Longitude must be between -180 and 180')
        
        self.latitude = latitude
        self.longitude = longitude
        self.time = time



    def get_zenith(self):
        """
        Get zenith coordinates
        """

        location = EarthLocation(lat=slef.latitude*u.deg, lon=self.longitude*u.deg)
        zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=AltAz(obstime=time, location=location))
        zenith_radec = zenith.transform_to('icrs')

        return zenith_radec.ra, zenith_radec.dec


test_object = User(latitude=48.85*u.deg, longitude=2.35*u.deg, time=Time.now())
print(test_object.latitude)

