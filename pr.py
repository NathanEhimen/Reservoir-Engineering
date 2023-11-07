import numpy as np 
import pandas as pd
import math

def calculate_AB(path_to_df):
    R = 8.31
    df = pd.read_csv(path_to_df, sep='\t') 
    a_mix = df['a_mix'].to_numpy()
    b_mix = df['b_mix'].to_numpy()
    P = df['P_Mpa'].to_numpy() * 10**6
    T = df['Temp'].to_numpy()

    A = a_mix * P / (R * T)**2 
    B = b_mix * P / R * T
    return np.vstack([A,B]).T


def PengRobinsonCubic(A, B,root=2):
        coeff = [A/A, B-1, A-2*B-3*B**2, B**3+B**2-A*B]
        a = coeff[0]
        b = coeff[1]
        c = coeff[2]
        d = coeff[3]
        
        p = (3 * a * c - b ** 2) / (3 * a ** 2)
        q = (2 * (b ** 3) - 9 * a * b * c + 27 * (a ** 2) * d) / (27 * a ** 3)
        discriminant = 18 * a * b * c * d - 4 * (b ** 3) * d + (b ** 2) * (c ** 2) - 4 * a * (c ** 3) - 27 * (a ** 2) * (d ** 2)
        
        if discriminant > 0: #there are 3 real roots
            return 2 * ((-p / 3) ** (1 / 2)) * np.cos((1 / 3) * np.arccos(((3 * q) / (2 * p)) * ((-3 / p) ** (1 / 2))) - 2 * np.pi * root / 3) - b / (3 * a)
        elif p > 0:
            return -2 * ((p / 3) ** (1 / 2)) * np.sinh((1 / 3) * np.arcsinh(((3 * q) / (2 * p)) * ((3 / p) ** (1 / 2)))) - b / (3 * a)
        else:
            return -2 * (np.abs(q) / q) * ((-p / 3) ** (1 / 2)) * np.cosh((1 / 3) * np.arccosh(((-3 * np.abs(q)) / (2 * p)) * ((-3 / p) ** (1 / 2)))) - b / (3 * a)


    # Index(['P_Mpa', 'ro_mol', 'Z', 'Temp', 'Mr', 'x_CO2', 'density',
    #    'mix_viscosity', 'a_mix', 'b_mix', 'p_red', 't_red'],