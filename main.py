from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord




class User(object):
    """
    User class
    """

    def __init__(self, ra, dec, time):
        """
        
        """
        self.ra = ra
        self.dec = dec
        self.time = Time(time)

    def get_zenith(self):
        """
        Get zenith coordinates
        """

        location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg)
        zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=AltAz(obstime=time, location=location))
        zenith_radec = zenith.transform_to('icrs')

        return zenith_radec.ra, zenith_radec.dec
