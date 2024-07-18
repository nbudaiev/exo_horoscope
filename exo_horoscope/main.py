from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import astropy.units as u
from geopy import geocoders  
from geopy.geocoders import Nominatim
import numpy as np

exoplanets_table = NasaExoplanetArchive.query_criteria(table="pscomppars", select="*")



class User(object):
    """
    User class
    """

    def __init__(self, user, citystate, year, month, day, hour, minute, second):
        """
        Args:
            user (str): User's name in the form 'User'
            citystate (str): City and State of birth in the form: 'City State' / 'City Country'
            year (int): Birthyear of User
            month (int): Birthmonth of User
            day (int): Birthday of User
            hour (int): Birthhour of User
            minute (int): Birthminute of User
            second (float): Birthsecond of User
        """

        if not isinstance(user, str):
            raise TypeError("User name must be a string.")

        if not isinstance(citystate, str):
            raise TypeError("City State / City Country must be a string.")

        if not isinstance(year, int):
            raise TypeError("Year must be an integer.")

        if year<=0:
            raise ValueError("Year must be a positive integer.")
        
        if not isinstance(month, int):
            raise TypeError("Month must be an integer.")
        if month<=0 or month>12:
            raise ValueError("Month must be an integer between 1 and 12.")
        
        if not isinstance(day, int):
            raise TypeError("Day must be an integer.")
        if day<=0 or day>31:
            raise ValueError("Day must be an integer between 1 and 31.")
        
        if not isinstance(hour, int):
            raise TypeError("Hour must be an integer.")
        if hour<0 or hour>23:
            raise ValueError("Hour must be an integer between 0 and 23.")
        
        if not isinstance(minute, int):
            raise TypeError("Minute must be an integer.")
        if minute<0 or minute>59:
            raise ValueError("Minute must be an integer between 0 and 59.")
        
        if not isinstance(second, (int, float)):
            raise TypeError("Second must be a float or an integer.")
        if second < 0 or second >= 60:
            raise ValueError("Second must be a float or an integer between 0 and 60.")
            raise TypeError("Second must be a float.")
        if second<0 or second>=60:
            raise ValueError("Second must be a float between 0 and 60.")

        

        self.user = user

        self.citystate = citystate

        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        time = Time(f'{self.year}-{self.month}-{self.day} {self.hour}:{self.minute}:{self.second}')
        self.time = time

        self.closest_object_nasa_table = self.get_closest_table()

    def get_closest_table(self):
        '''
        Get table of closest object

        This method finds the Nasa Exoplanet Archive table of the object which transits nearest birth zenith of the user.

        Args:
            citystate (str): City and State of birth in the form: 'City State'
            year (int): Birthyear of User
            month (int): Birthmonth of User
            day (int): Birthday of User
            hour (int): Birthhour of User
            minute (int): Birthminute of User
            second (float): Birthsecond of User

        Returns:
            astropy.table.table.QTable: table of closest object to birth zenith
        '''
            

        geolocator = Nominatim(user_agent='moeur')
        location = geolocator.geocode(self.citystate)
        self.birth_lat, self.birth_lon = location[1][0], location[1][1]
        table = exoplanets_table
        coords = self.get_zenith()
        stars_coords = SkyCoord(table['ra'], table['dec'], unit=(u.deg, u.deg))
        distances = coords.separation(stars_coords)
        closest_index = distances.argmin()
        closest_object = table[closest_index]
        closest_table = NasaExoplanetArchive.query_object(closest_object['pl_name'])
        return closest_table
        
    def get_zenith(self):
        """
        Compute birth zenith

        This method takes latitude and longitude coordinates of the user's birth city and time of birth and calculates the celestial coordinates of the zenith.
        
        Returns:
            astropy.coordinates.sky_coordinate.SkyCoord: celestial coordinates of the zenith.
        """
        location = EarthLocation(lat=self.birth_lat, lon=self.birth_lon)
        zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=AltAz(obstime=self.time, location=location))
        zenith_radec = zenith.transform_to('icrs')

        return zenith_radec



    def map_eccentricity_to_trait(self):
        '''
        Map orbital eccentricity to personality trait

        This method assigns a personality trait to the user based on the value of their birth exoplanet's orbital eccentricity.

        Returns:
            str: the personality trait
        '''
        if self.eccentricity == np.nan:
            return ""
        if self.eccentricity == 0:
            return "are perfectly stable"
        elif 0 < self.eccentricity < 0.3:
            return "prefer stability"
        elif 0.3 <= self.eccentricity < 0.6:
            return "are balanced"
        elif 0.6 <= self.eccentricity < 0.9:
            return "prefer excitement"
        else:
            return "embrace change"

    def map_semimajor_axis_to_trait(self):
        '''
        Map orbital semimajor axis to personality trait

        This method assigns a personality trait to the user based on the value of their birth exoplanet's orbital semimajor axis.

        Returns:
            str: the personality trait
        '''
        if self.semimajor_axis == np.nan:
            return ""
        if self.semimajor_axis < 0.1:
            return "are 'close to the action' and constantly influenced by your star's energy, suggesting a very outgoing and active nature"
        elif 0.1 <= self.semimajor_axis < 1:
            return "are still within a region of significant stellar influence, indicating a generally social and engaging character"
        elif 1 <= self.semimajor_axis < 5:
            return "strike a balance between the inner and outer regions, reflecting a well-rounded personality that is equally comfortable in social situations and solitude"
        elif 5 <= self.semimajor_axis < 30:
            return "are farther from the star, implying a more reserved and introspective nature, preferring less direct interaction"
        else:
            return "are on the outskirts, indicating a highly introspective and solitary disposition, thriving in their own space away from the hustle and bustle"

    def map_orbital_period_to_trait(self):
        """
        Map exoplanet system's orbital period to a personality trait (thinking style).

        Returns: 
            str: The personality trait text message.
        """
        if self.period == np.nan:
            return ""
        if self.period < 10:
            return "rapid orbits suggest a fast-paced and reactive thinking style"
        elif 10 <= self.period < 100:
            return "orbits allow for rapid changes and adaptation, indicating an active and adaptive thinking style"
        elif 100 <= self.period < 365:
            return "orbital periods tend to experience balanced conditions, suggesting a balanced and analytical thinking style"
        elif 365 <= self.period < 3650:
            return "planets take longer to orbit their stars, implying a more deliberate and thoughtful approach"
        else:
            return "very long orbital periods embody a reflective and contemplative thinking style"

    def map_stellar_mass_to_trait(self):
        """
        Map exoplanet system's stellar mass to a personality trait.

        Returns:
            str: The personality trait based on the stellar mass.
        """

        if self.stellar_mass == np.nan:
            return ""
        if self.stellar_mass < 0.5:
            return "stable and enduring"
        elif 0.5 <= self.stellar_mass < 1.5:
            return "balanced and nurturing"
        elif 1.5 <= self.stellar_mass < 3:
            return "dynamic and charismatic"
        else:
            return "intense and transformative"
        
    def get_horoscope(self):
        """
        User class method to get the User's horoscope based on User's attributes.

        Returns:
            str: The horoscope message for the User.
        """
        self.planet = self.closest_object_nasa_table['pl_name'][0]
        self.star = self.closest_object_nasa_table['hostname'][0]
        self.eccentricity = np.nanmean(self.closest_object_nasa_table["pl_orbeccen"])
        self.semimajor_axis = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_orbsmax"].value))
        self.period = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_orbper"].value))
        self.stellar_mass = np.nanmean(np.asarray(self.closest_object_nasa_table["st_mass"].value))
        #self.planet_mass = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_bmassj"].value))
        #self.planet_radius = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_radj"].value))
        #self.planet_density = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_dens"].value))
        #self.planet_temp = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_eqt"].value))
        #self.stellar_magnitude = np.nanmean(np.asarray(self.closest_object_nasa_table["st_optmag"].value))
        #self.stellar_radius = np.nanmean(np.asarray(self.closest_object_nasa_table["st_rad"].value))
        #self.stellar_temp = np.nanmean(np.asarray(self.closest_object_nasa_table["st_teff"].value))


        eccentricity_trait = self.map_eccentricity_to_trait()
        axis_trait = self.map_semimajor_axis_to_trait()
        period_trait = self.map_orbital_period_to_trait()
        stellar_mass_trait = self.map_stellar_mass_to_trait()

        message = (f"{self.user}, your birth exoplanet is {self.planet} orbiting star {self.star}. "
                f"Based on an eccentricity of {self.eccentricity:.2f}, you {self.map_eccentricity_to_trait()}. "
                f"With an orbit semi-major axis of {self.semimajor_axis:.2f} AU, you {self.map_semimajor_axis_to_trait()}. "
                f"With a birth exoplanet period of {self.period:.2f} days, these {self.map_orbital_period_to_trait()}, "
                f"and with a stellar mass of {self.stellar_mass:.2f} solar masses, you are {self.map_stellar_mass_to_trait()}.")
        return message

