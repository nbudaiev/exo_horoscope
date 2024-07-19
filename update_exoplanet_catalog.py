from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

# load the current catalog of confirmed exoplanets
exoplanets_table = NasaExoplanetArchive.query_criteria(table="pscomppars", select="*")

exoplanets_table_selected_columns = exoplanets_table['pl_name', 'hostname','pl_orbeccen', 'pl_orbsmax','pl_orbper','st_mass',
                                                    'pl_bmassj', 'pl_radj', 'pl_dens', 'pl_eqt', 'st_rad', 'st_teff']

exoplanets_table_selected_columns.write("confirmed_exoplanets_table.ecsv", format="ascii.ecsv", overwrite=True)