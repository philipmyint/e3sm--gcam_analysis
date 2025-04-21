import numpy as np

# Conversion factors.
km_TO_m = 1e3
km2_TO_m2 = km_TO_m*km_TO_m
days_TO_s = 86400
kg_TO_g = 1e3
Pg_TO_kg = 1e12
Pg_TO_g = Pg_TO_kg*kg_TO_g
m_TO_mm = 1e3
mole_fraction_TO_ppm = 1e6    

# Number of days and seconds in each month (no leap years).
NUM_DAYS_IN_MONTHS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
NUM_SECONDS_IN_MONTHS = days_TO_s*NUM_DAYS_IN_MONTHS

# Molar mass of carbon (g/mol).
MM_C = 12.011

# Molar mass of CO2 (g/mol).
MM_CO2 = 44.009

# Molar mass of CH4 (g/mol).
MM_CH4 = 16.043

# Molar mass of atmosphere avg (g/mol) (dry).
MM_ATM = 28.965

# Surface area of Earth in m^2. Based on authalic radius, which is spherical radius that gives same 
# surface area as the reference ellipsoid, which is wgs84.
SURF_AREA = 4 * (4.0 * np.atan(1.0)) * 6371007 * 6371007