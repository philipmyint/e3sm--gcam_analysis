import numpy as np

""" Unit conversion factors. """
km_TO_m = 1e3
km2_TO_m2 = km_TO_m*km_TO_m
days_TO_s = 86400
days_TO_hours = 24
hours_TO_mins = 60
years_TO_days = 365
years_TO_s = years_TO_days*days_TO_s
years_TO_months = 12
years_TO_hours = years_TO_days*days_TO_hours
years_TO_mins = years_TO_hours*hours_TO_mins
kg_TO_g = 1e3
Pg_TO_kg = 1e12
Pg_TO_g = Pg_TO_kg*kg_TO_g
m_TO_mm = 1e3
mole_fraction_TO_ppm = 1e6    

""" Tiny constant used to avoid divide-by-zero errors. """
EPSILON = 1.0e-14

""" Number of days and seconds in each month (no leap years). """
NUM_DAYS_IN_MONTHS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
NUM_SECONDS_IN_MONTHS = days_TO_s*NUM_DAYS_IN_MONTHS

""" Molar mass of carbon (g/mol). """
MM_C = 12.011

""" Molar mass of CO2 (g/mol). """
MM_CO2 = 44.009

""" Molar mass of CH4 (g/mol). """
MM_CH4 = 16.043

""" Molar mass of atmosphere avg (g/mol) (dry). """
MM_ATM = 28.965

""" Surface area of Earth in m^2. Based on authalic radius, which is a spherical radius that gives same surface area 
as the reference ellipsoid, which is wgs84. """
SURF_AREA = 4 * (4.0 * np.atan(1.0)) * 6371007 * 6371007

""" Dictionary between month number and name. """
MONTH_NUM_TO_NAME = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June',
                     7: 'July', 8: 'August', 9: 'September', 10: 'October', 11: 'November', 12: 'December'}