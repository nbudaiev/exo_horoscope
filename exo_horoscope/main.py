from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import astropy.units as u
from geopy.geocoders import Nominatim
import numpy as np
from astropy.io import ascii
import warnings
import importlib.resources

import os

#file_path = "confirmed_exoplantes_table.ecsv"

with importlib.resources.path('exo_horoscope', 'main.py') as package_root_path:
    package_root = package_root_path.parent

catalog_path = os.path.join(package_root, 'confirmed_exoplanets_table.ecsv')

if not os.path.exists(catalog_path):
   from exo_horoscope import update_exoplanet_catalog

with importlib.resources.path('exo_horoscope', 'confirmed_exoplanets_table.ecsv') as catalog_path:
    exoplanets_table = ascii.read(catalog_path)



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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            date_and_time = Time(f'{self.year}-{self.month}-{self.day} {self.hour}:{self.minute}:{self.second}')
        self.time = date_and_time
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.closest_object_nasa_table = self.get_closest_table()

        self.planet = self.closest_object_nasa_table['pl_name']
        self.star = self.closest_object_nasa_table['hostname']



    def get_closest_table(self):
        """
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
        """
            
        geolocator = Nominatim(user_agent='moeur')
        location = geolocator.geocode(self.citystate)
        self.birth_lat, self.birth_lon = location[1][0], location[1][1]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coords = self.get_zenith()
            stars_coords = SkyCoord(exoplanets_table['ra'], exoplanets_table['dec'], unit=(u.deg, u.deg))
        distances = coords.separation(stars_coords)
        closest_index = distances.argmin()
        closest_table = exoplanets_table[closest_index]
        return closest_table
        
    def get_zenith(self):
        """
        Compute birth zenith

        This method takes latitude and longitude coordinates of the user's birth city and time of birth and calculates the celestial coordinates of the zenith.
        
        Returns:
            astropy.coordinates.sky_coordinate.SkyCoord: celestial coordinates of the zenith.
        """
        location = EarthLocation(lat=self.birth_lat, lon=self.birth_lon)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=AltAz(obstime=self.time, location=location))
            zenith_radec = zenith.transform_to('icrs')
        return zenith_radec



    def map_eccentricity_to_trait(self):
        """
        Map orbital eccentricity to personality trait

        This method assigns a personality trait to the user based on the value of their birth exoplanet's orbital eccentricity.

        Returns:
            str: the personality trait
        """
        if np.isnan(self.eccentricity):
            return ""

        base_traits = [
            "are perfectly stable, thriving in predictable and consistent environments",
            "have a very low tendency for change, enjoying a calm and steady lifestyle",
            "exhibit a slight preference for stability, appreciating consistency",
            "favor routine, enjoying the security of familiar patterns",
            "are mostly stable, finding comfort in regularity",
            "enjoy a generally stable environment, feeling at ease with the familiar",
            "have a low eccentricity, preferring a predictable and steady life",
            "lean towards stability, valuing consistency and routine",
            "enjoy a slightly varied life, with a preference for stability",
            "are stable but open to small changes, appreciating minor variations",
            "prefer stability but enjoy occasional changes, valuing a balance",
            "are generally stable, with a slight openness to new experiences",
            "enjoy a balanced life, appreciating both stability and change",
            "are balanced, finding comfort in a mix of stability and variety",
            "appreciate a balanced approach, valuing both routine and novelty",
            "enjoy equilibrium, equally comfortable with consistency and change",
            "are balanced, thriving in environments with both stability and excitement",
            "value a balanced life, enjoying a mix of familiar and new experiences",
            "are balanced, appreciating the interplay of stability and novelty",
            "prefer excitement but value stability, enjoying a dynamic life",
            "enjoy a mix of excitement and stability, thriving in varied environments",
            "lean towards excitement, appreciating dynamic and stimulating experiences",
            "prefer excitement, valuing variety and change",
            "are drawn to excitement, thriving in lively and dynamic settings",
            "enjoy excitement, preferring stimulating and varied experiences",
            "prefer a dynamic life, valuing excitement and change",
            "are highly dynamic, thriving in exciting and varied environments",
            "embrace change, enjoying a life full of variety and excitement",
            "are very dynamic, appreciating continuous change and excitement",
            "embrace change, thriving in unpredictable and dynamic environments"
        ]
        
        # Define the number of bins and the upper bound for eccentricity
        num_bins = len(base_traits)
        upper_bound = 1.0  # Eccentricity ranges from 0 to 1

        # Determine the bin width
        bin_width = upper_bound / num_bins
        
        # Determine the bin index for the eccentricity
        bin_index = int(self.eccentricity / bin_width)
        
        # Ensure the index is within the bounds of base_traits
        bin_index = min(bin_index, num_bins - 1)
        
        return base_traits[bin_index]

    def map_semimajor_axis_to_trait(self):
        """
        Map orbital semimajor axis to personality trait

        This method assigns a personality trait to the user based on the value of their birth exoplanet's orbital semimajor axis.

        Returns:
            str: the personality trait
        """
        if np.isnan(self.semimajor_axis):
            return ""

        base_traits = [
            "are extremely close to the action, feeling the intensity of your star's energy, suggesting a very outgoing and hyperactive nature",
            "are very near your star, constantly in the spotlight, indicating a lively and energetic personality",
            "experience strong stellar influences, reflecting a dynamic and enthusiastic character",
            "are close enough to feel the heat, suggesting a vibrant and extroverted nature",
            "are quite near your star, basking in its energy, implying a spirited and sociable disposition",
            "feel the star's warmth intensely, reflecting an animated and vivacious personality",
            "are within close range, constantly energized by your star, suggesting a bubbly and proactive nature",
            "are near the stellar core, indicating a dynamic and gregarious character",
            "receive strong stellar input, implying an outgoing and lively personality",
            "are influenced heavily by your star, reflecting an exuberant and active disposition",
            "are still within a region of significant stellar influence, indicating a generally social and engaging character",
            "are within the star's active zone, suggesting a sociable and energetic personality",
            "are close enough for strong stellar effects, implying a friendly and dynamic nature",
            "are in the star's warm embrace, reflecting a warm and approachable personality",
            "experience significant stellar interaction, suggesting a vibrant and engaging character",
            "are in a region of notable stellar influence, indicating a lively and amiable nature",
            "are within the star's reach, suggesting a sociable and spirited personality",
            "are moderately influenced by the star, reflecting a balanced and outgoing nature",
            "are within comfortable stellar proximity, indicating a well-rounded and friendly character",
            "are in a balanced zone, implying a sociable and dynamic personality",
            "strike a balance between the inner and outer regions, reflecting a well-rounded personality that is equally comfortable in social situations and solitude",
            "are in a region of moderate stellar influence, suggesting a balanced and thoughtful nature",
            "are at a moderate distance from the star, indicating a balanced and introspective personality",
            "are slightly removed from the star, reflecting a balanced and contemplative character",
            "are somewhat distanced from the star, suggesting a reserved and thoughtful nature",
            "are farther from the star, implying a more reserved and introspective nature, preferring less direct interaction",
            "are in a region of lesser stellar influence, indicating a reserved and reflective personality",
            "are at a significant distance, reflecting a thoughtful and solitary nature",
            "are quite removed from the star, suggesting a contemplative and independent character",
            "are on the outskirts, indicating a highly introspective and solitary disposition, thriving in their own space away from the hustle and bustle"
        ]
        
        # Define the number of bins and the upper bound for semimajor axis
        num_bins = len(base_traits)
        upper_bound = 30.0  # Assuming 30 AU as the upper bound for semimajor axis

        # Determine the bin width
        bin_width = upper_bound / num_bins
        
        # Determine the bin index for the semimajor axis
        bin_index = int(self.semimajor_axis / bin_width)
        
        # Ensure the index is within the bounds of base_traits
        bin_index = min(bin_index, num_bins - 1)
        
        return base_traits[bin_index]

    def map_orbital_period_to_trait(self):
        """
        Map exoplanet system's orbital period to a personality trait (thinking style).

        Returns:
            str: The personality trait text message.
        """
        if np.isnan(self.period):
            return ""
        
        base_traits = [
            "rapid orbits suggest a fast-paced and reactive thinking style",
            "quick orbits imply a highly adaptable and spontaneous nature",
            "short orbits point to an energetic and proactive approach",
            "swift orbits indicate a keen sense of urgency and decisiveness",
            "brisk orbits reflect a sharp and agile mind",
            "speedy orbits suggest a dynamic and flexible thinking style",
            "fast orbits imply a quick-witted and responsive nature",
            "zippy orbits point to a lively and alert mentality",
            "nippy orbits indicate a nimble and resourceful mind",
            "hasty orbits reflect a prompt and efficient thinking style",
            "moderate orbits allow for rapid changes and adaptation, indicating an active and adaptive thinking style",
            "steady orbits suggest a balanced and rational approach",
            "intermediate orbits point to a thoughtful and considerate nature",
            "consistent orbits imply a stable and reliable thinking style",
            "reliable orbits reflect a methodical and logical approach",
            "balanced orbits tend to experience balanced conditions, suggesting a balanced and analytical thinking style",
            "measured orbits imply a careful and precise nature",
            "calm orbits point to a composed and reflective mindset",
            "harmonious orbits suggest an even-tempered and insightful approach",
            "stable orbits reflect a judicious and perceptive thinking style",
            "long orbits take longer to orbit their stars, implying a more deliberate and thoughtful approach",
            "extended orbits suggest a contemplative and introspective nature",
            "prolonged orbits point to a patient and reflective mindset",
            "lengthy orbits imply a thorough and meticulous approach",
            "drawn-out orbits reflect a cautious and prudent nature",
            "vast orbits embody a reflective and contemplative thinking style",
            "sprawling orbits suggest a far-reaching and philosophical mindset",
            "expansive orbits point to a profound and meditative nature",
            "sweeping orbits imply a deep and contemplative approach",
            "broad orbits reflect an insightful and discerning thinking style"
        ]
        
        # Define the number of bins and the upper bound for orbital periods
        num_bins = len(base_traits)
        upper_bound = 3650.0  # Assuming 3650 days as the upper bound for orbital periods

        # Determine the bin width
        bin_width = upper_bound / num_bins
        
        # Determine the bin index for the orbital period
        bin_index = int(self.period / bin_width)
        
        # Ensure the index is within the bounds of base_traits
        bin_index = min(bin_index, num_bins - 1)
        
        return base_traits[bin_index]

    def map_stellar_mass_to_trait(self):
        """
        Map exoplanet system's stellar mass to a personality trait.

        Returns:
            str: The personality trait based on the stellar mass.
        """
        if np.isnan(self.stellar_mass):
            return ""

        base_traits = [
            "stable and enduring", "patient and reliable", "balanced and nurturing",
            "adaptable and flexible", "thoughtful and observant", "resilient and steady",
            "meticulous and precise", "practical and grounded", "loyal and trustworthy",
            "persistent and determined", "independent and self-sufficient", "resourceful and clever",
            "curious and inquisitive", "imaginative and dreamy", "empathetic and compassionate",
            "charming and persuasive", "versatile and multifaceted", "playful and spontaneous",
            "innovative and creative", "dynamic and charismatic", "confident and bold",
            "ambitious and driven", "resolute and steadfast", "calm and composed",
            "energetic and lively", "visionary and forward-thinking", "wise and discerning",
            "courageous and daring", "intense and transformative", "powerful and dominant"
        ]
        
        # Define the number of bins and the upper bound for stellar mass
        num_bins = len(base_traits)
        upper_bound = 3.0  # Assuming 3.0 as the upper bound for stellar mass

        # Determine the bin width
        bin_width = upper_bound / num_bins
        
        # Determine the bin index for the stellar mass
        bin_index = int(self.stellar_mass / bin_width)
        
        # Ensure the index is within the bounds of base_traits
        bin_index = min(bin_index, num_bins - 1)
        
        return base_traits[bin_index]

        
    def get_horoscope(self):
        """
        User class method to get the User's horoscope based on User's attributes.

        Returns:
            str: The horoscope message for the User.
        """
        self.eccentricity = np.nanmean(self.closest_object_nasa_table["pl_orbeccen"])
        self.semimajor_axis = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_orbsmax"].value))
        self.period = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_orbper"].value))
        self.stellar_mass = np.nanmean(np.asarray(self.closest_object_nasa_table["st_mass"].value))
        #self.planet_mass = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_bmassj"].value))
        #self.planet_radius = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_radj"].value))
        #self.planet_density = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_dens"].value))
        #self.planet_temp = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_eqt"].value))
        #self.stellar_magnitude = np.nanmean(np.asarray(self.closest_object_nasa_table["st_optmag"].value)) # I think this one is wrong
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


    def map_radius_to_life_suggestion(self):
        """
        Map planet radius to life suggestion.

        This method assigns a life suggestion to the user based on the value of their birth exoplanet's radius.

        Returns:
            str: the life suggestion
        """
        if np.isnan(self.radius):
            return ""

        base_suggestions = [
            "focus on the little things; small steps can lead to big achievements",
            "pay attention to detail; even small efforts can have significant impact",
            "embrace small changes; they can lead to greater success",
            "take gradual steps towards your goals; patience will be rewarding",
            "be attentive to minor aspects; they often hold the key to success",
            "celebrate small victories; they pave the way for bigger accomplishments",
            "appreciate the small things; they contribute to overall growth",
            "find joy in small achievements; they build the path to larger goals",
            "take small, consistent actions; they lead to significant outcomes",
            "value incremental progress; it accumulates to substantial results",
            "balance ambition with contentment; focus on steady progress",
            "set realistic goals; steady progress leads to achievement",
            "maintain equilibrium between ambition and satisfaction",
            "find contentment in your journey; progress comes with patience",
            "cultivate a sense of balance; aim for steady and sustainable growth",
            "focus on maintaining stability; small successes build up over time",
            "find harmony between your aspirations and achievements",
            "embrace a balanced approach; gradual progress is valuable",
            "be mindful of your progress; find satisfaction in steady growth",
            "find a middle ground; balance your drive with appreciation for progress",
            "be bold and take on challenges head-on; step out of your comfort zone",
            "pursue opportunities with courage; face challenges directly",
            "take initiative; don't hesitate to embrace new challenges",
            "be proactive in pursuing your goals; confront challenges with confidence",
            "tackle obstacles with determination; take bold steps forward",
            "embrace challenges; they are opportunities for growth and learning",
            "pursue ambitious goals; your boldness will lead to growth",
            "take on new ventures with confidence; challenges are stepping stones",
            "be assertive in your ambitions; tackle obstacles with enthusiasm",
            "aim high and don't be afraid to dream big; pursue your grand visions",
            "aspire to great heights; let your ambitions drive you forward",
            "dream without limits; your aspirations are the key to success",
            "set lofty goals; your determination will help you achieve them",
            "reach for the stars; embrace ambitious pursuits and grand visions",
            "strive for greatness; let your dreams inspire you to new heights"
        ]
        
        num_bins = len(base_suggestions)
        upper_bound = 10.0  # Assuming 10 as the upper bound for radius

        bin_width = upper_bound / num_bins
        bin_index = int(self.radius / bin_width)
        bin_index = min(bin_index, num_bins - 1)

        return base_suggestions[bin_index]

    def map_magnitude_to_life_suggestion(self):
        """
        Map planet magnitude to life suggestion.

        This method assigns a life suggestion to the user based on the value of their birth exoplanet's magnitude.

        Returns:
            str: the life suggestion
        """
        if np.isnan(self.magnitude):
            return ""

        base_suggestions = [
            "embrace your bright and positive nature",
            "let your inner light shine brightly; positivity is your strength",
            "radiate positivity; your brightness inspires others",
            "share your optimistic outlook; it illuminates those around you",
            "maintain your radiant energy; it enhances your impact",
            "celebrate your natural brightness; it brings light to others",
            "be a beacon of positivity; your brightness is inspiring",
            "find ways to shine even in the dark moments; your resilience is key",
            "bring light to challenging situations; your optimism is powerful",
            "continue to shine through adversity; your brightness leads the way",
            "spread your light in difficult times; it has a profound effect",
            "help others see the light; your positivity can guide them",
            "be a guiding light for others around you; your presence is illuminating",
            "lead by example; your brightness provides direction and hope",
            "use your light to help others navigate through their challenges",
            "offer guidance with your bright perspective; it lights the path for others",
            "inspire those around you with your bright presence; it offers hope",
            "be a source of illumination; your brightness has a meaningful impact",
            "shine brightly; your presence brings clarity and inspiration",
            "be the light in the lives of others; your brightness is a gift",
            "navigate challenges with your inner light; it guides you through",
            "find strength in your ability to shine; it lights up even the darkest moments",
            "offer guidance and support with your illuminating presence",
            "be a source of warmth and light; your positivity is transformative",
            "continue to be a beacon of hope; your brightness encourages growth",
            "let your light shine unwaveringly; it creates a path for others",
            "radiate warmth and positivity; your light has a profound effect",
            "be a steady source of illumination; your presence inspires and guides",
            "embrace your role as a guiding light; your brightness fosters growth"
        ]
        
        num_bins = len(base_suggestions)
        upper_bound = 20.0  # Assuming 20 as the upper bound for magnitude

        bin_width = upper_bound / num_bins
        bin_index = int(self.magnitude / bin_width)
        bin_index = min(bin_index, num_bins - 1)

        return base_suggestions[bin_index]

    def map_density_to_life_suggestion(self):
        """
        Map planet density to life suggestion.

        This method assigns a life suggestion to the user based on the value of their birth exoplanet's density.

        Returns:
            str: the life suggestion
        """
        if np.isnan(self.density):
            return ""

        base_suggestions = [
            "keep a light-hearted and flexible approach to life",
            "maintain a flexible attitude; adapt easily to changes",
            "embrace a light-hearted approach; it brings joy and ease",
            "find joy in being adaptable; flexibility is your strength",
            "stay open to changes; flexibility will enhance your experiences",
            "approach life with a light touch; adaptability is key",
            "keep a carefree attitude; it helps you navigate through challenges",
            "balance flexibility with stability; it leads to a fulfilling life",
            "embrace adaptability; it makes life's journey more enjoyable",
            "be light-hearted and adaptable; it helps you handle various situations",
            "balance your seriousness with moments of joy; it enriches your life",
            "find joy amidst your responsibilities; it makes your approach more enjoyable",
            "incorporate fun into your serious tasks; it brings balance",
            "blend seriousness with light-hearted moments; it creates harmony",
            "maintain a balance between work and play; it enhances your well-being",
            "find joy in the midst of your responsibilities; it brings balance",
            "combine diligence with moments of relaxation; it improves your outlook",
            "find equilibrium between seriousness and enjoyment; it leads to satisfaction",
            "balance your approach to life; mix seriousness with playfulness",
            "stay grounded and practical in your decisions; it provides stability",
            "be practical and realistic; it ensures steady progress",
            "maintain a grounded approach; it helps you make effective decisions",
            "focus on practical solutions; it leads to reliable outcomes",
            "remain grounded in your approach; it supports long-term success",
            "be practical and pragmatic; it keeps you focused on achievable goals",
            "be resilient and unyielding in the face of challenges; it ensures progress",
            "develop resilience; it helps you overcome obstacles effectively",
            "approach challenges with determination; resilience will see you through",
            "remain steadfast in adversity; your resilience is a key strength",
            "be unwavering in your efforts; resilience will lead to success",
            "build your resilience; it enables you to handle challenges effectively",
            "maintain a strong resolve; it supports your ability to overcome difficulties",
            "embrace your resilient nature; it empowers you to face challenges head-on",
            "stay firm in your approach; resilience ensures progress despite obstacles"
        ]
        
        num_bins = len(base_suggestions)
        upper_bound = 10.0  # Assuming 10 as the upper bound for density

        bin_width = upper_bound / num_bins
        bin_index = int(self.density / bin_width)
        bin_index = min(bin_index, num_bins - 1)

        return base_suggestions[bin_index]

    def get_life_suggestions(self):
        """
        User class method to get the User's life suggestions based on exoplanet's orbital properties.

        Returns:
            str: The life suggestions message for the User.
        """
        self.radius = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_radj"].value))
        self.magnitude = np.nanmean(np.asarray(self.closest_object_nasa_table["sy_gaiamag"].value))
        self.density = np.nanmean(np.asarray(self.closest_object_nasa_table["pl_dens"].value))
        radius_suggestion = self.map_radius_to_life_suggestion()
        magnitude_suggestion = self.map_magnitude_to_life_suggestion()
        density_suggestion = self.map_density_to_life_suggestion()

        message = (f"{self.user}, your birth exoplanet is {self.planet} orbiting star {self.star}. "
                f"Based on a radius of {self.radius:.2f} Jupiter radii, {radius_suggestion}. "
                f"With a magnitude of {self.magnitude:.2f}, {magnitude_suggestion}. "
                f"And with a density of {self.density:.2f} g/cm³, {density_suggestion}.")
        return message
    
    def get_lucky_numbers(self):
        """
        Generate lucky numbers based on the first two letters of the exoplanet and user names.

        Returns:
            str: A message with the lucky numbers and their corresponding adjectives.
        """
        # Dictionary to map letters to their positions in the alphabet
        letter_to_number = {chr(i + 96): i for i in range(1, 27)}

        # Function to get number from letter
        def letter_number(letter):
            return letter_to_number.get(letter.lower(), 0)

        # Function to get an adjective based on a number
        def number_to_adjective(number):
            adjectives = [
                "amazing", "brave", "creative", "dynamic", "elegant", "fearless", "graceful", "honest",
                "intelligent", "joyful", "kind", "lively", "mighty", "noble", "optimistic", "passionate",
                "quick", "radiant", "strong", "trustworthy", "unique", "vibrant", "wise", "youthful", "zealous"
            ]
            return adjectives[number % len(adjectives)]

        # Get the first two letters of the exoplanet and user names
        planet_letters = self.planet[:2].lower()
        user_letters = self.user[:2].lower()

        # Generate the lucky numbers
        lucky_numbers = [letter_number(letter) for letter in planet_letters + user_letters]

        # Generate the message with lucky numbers and adjectives
        lucky_numbers_message = ", ".join([f"{num} ({number_to_adjective(num)})" for num in lucky_numbers])

        message = f"Lucky numbers: {lucky_numbers_message}."
        return message