from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

# load the current catalog of confirmed exoplanets
exoplanets_table = NasaExoplanetArchive.query_criteria(table="pscomppars", select="*")

exoplanets_table.write("confirmed_exoplanets_table_test.ecsv", format="ascii.ecsv", overwrite=True)