#from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import astropy.units as u
from geopy import geocoders  
from geopy.geocoders import Nominatim
from datetime import datetime

class User(object):
    """
    User class
    """

    def __init__(self, user, citystate, year, month, day, hour, minute, second):
        """
        Args:
        name (str): Name of User
        citystate (str): City and State of birth in the form: 'City State'
        year (int): Birthyear of User
        month (int): Birthmonth of User
        day (int): Birthday of User
        hour (int): Birthhour of User
        minute (int): Birthminute of User
        second (int): Birthsecond of User
        """
        self.user = user
        self.citystate = citystate
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.closest_object_nasa_table = self.get_closest_table()

    def get_closest_table(self):
            
        geolocator = Nominatim(user_agent='moeur')
        location = geolocator.geocode(self.citystate)
        date_time = datetime(self.year, self.month, self.day, self.hour, self.minute, self.second)
        birth_ra, birth_dec = location[1][0], location[1][1]
        zen_ra, zen_dec = zenith_ra_dec(birth_ra, birth_dec, date_time)
        table = NasaExoplanetArchive.query_region(table="ps", coordinates=SkyCoord(zen_ra * u.deg, zen_dec * u.deg), radius=4 * u.deg)
        coords = SkyCoord(zen_ra * u.deg, zen_dec * u.deg)
        stars_coords = SkyCoord(table['ra'], table['dec'], unit=(u.deg, u.deg))
        distances = coords.separation(stars_coords)
        closest_index = distances.argmin()
        closest_object = table[closest_index]
        closest_table = NasaExoplanetArchive.query_object(closest_object['pl_name'])
        return closest_table
        
        
    def get_horoscope(self):
        self.closest_object_nasa_table['pl_name'][0]
        star = self.closest_object_nasa_table['hostname'][0]
        eccentricity = np.nanmean(self.closest_object_nasa_table["pl_orbeccen"])
        semi_major_axis = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_orbsmax"].value))
        period = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_orbper"].value))
        stellar_mass = np.nanmean(np.asarray(self.closest_object_nasa_table["st_mass"].value))
        eccentricity_trait = map_eccentricity_to_trait(eccentricity)
        axis_trait = map_semimajor_axis_to_trait(semi_major_axis)
        period_trait = map_orbital_period_to_trait(period)
        stellar_mass_trait = map_stellar_mass_to_trait(stellar_mass)

        message = (f"{self.user}, your birth exoplanet is {name} orbiting star {star}. "
                f"Based on an eccentricity of {eccentricity:.2f}, you {eccentricity_trait}. "
                f"With an orbit semi-major axis of {semi_major_axis:.2f} AU, you {axis_trait}. "
                f"With a birth exoplanet period of {period:.2f} days, these {period_trait}, "
                f"and with a stellar mass of {stellar_mass:.2f} solar masses, you are {stellar_mass_trait}.")
        return message

#chelsea = User("Chelsea", "Redlands California", 1991, 4, 27, 12, 4, 0).get_horoscope()
        #geolocator = Nominatim(user_agent='moeur')
        #location = geolocator.geocode(citystate)
        #birth_lat, birth_log = location[1][0], location[1][1]
        #return birth_lat, birth_log
        #if latitude.unit != u.deg:
        #    raise ValueError('Latitude must be in degrees')
        #if longitude.unit != u.deg:
        #    raise ValueError('Longitude must be in degrees')

        #if latitude.value < -90 or latitude.value > 90:
        #    raise ValueError('Latitude must be between -90 and 90')
        #if longitude.value < -180 or longitude.value > 180:
        #    raise ValueError('Longitude must be between -180 and 180')
        
        #self.latitude = latitude
        #self.longitude = longitude
        
        #self.time = time
        # implement the time using year, month, day,... version (similar to below.) Change the __init__ parameters
        #self.time = Time(datetime(year, month, day, hours, minute, second))
    
    
        #self.ra, self.dec = self.get_zenith()




    # def get_zenith(self):
    #    """
    #    Get zenith coordinates
    #    """

    #    location = EarthLocation(lat=self.latitude, lon=self.longitude)
    #    zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=AltAz(obstime=self.time, location=location))
    #    zenith_radec = zenith.transform_to('icrs')



        return zenith_radec.ra, zenith_radec.dec


#test_object = User(48.85*u.deg, 2.35*u.deg, Time.now())
#print(test_object.latitude)
