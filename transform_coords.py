import numpy as np
import pandas as pd

def cos(angle):
    return np.cos(np.deg2rad(float(angle)))

def sin(angle):
    return np.sin(np.deg2rad(float(angle)))

def get_T():
    ra_ngp = 192.7278
    dec_ngp = 26.8630
    theta = 122.9280
    T1 = np.array([[cos(theta), sin(theta), 0],
          [sin(theta), -cos(theta), 0],
          [0, 0, 1]])
    T2 = np.array([[-sin(dec_ngp), 0, cos(dec_ngp)],
          [0, -1, 0],
          [cos(dec_ngp), 0, sin(dec_ngp)]])
    T3 = np.array([[cos(ra_ngp), sin(ra_ngp), 0],
          [sin(ra_ngp), -cos(ra_ngp), 0],
          [0, 0, 1]])
    T = T1 @ T2 @ T3
    return T

def transform_coords(ra, dec, T):
    M1 = np.array([[cos(dec) * cos(ra)],
                   [cos(dec) * sin(ra)],
                   [sin(dec)]])
    M = T @ M1

    b = np.rad2deg(np.arcsin(M[2]))
    sin_l = M[0]/cos(b)
    l = np.rad2deg(np.arctan(M[1]/M[0]))
    if sin_l < 0:
        l += 180
    if l < 0:
        l = 360 + l
    return (l, b)


def get_velocity(ra, dec, R_v, mu_ra_cosdec, mu_dec, R, T):
    k = 4.740470463533

    A = (
            np.array([[cos(ra), sin(ra), 0],
                  [sin(ra), -cos(ra), 0],
                  [0, 0, -1]])
            @
            np.array([[cos(dec), 0, -sin(dec)],
                      [0, -1, 0],
                      [-sin(dec), 0, -cos(dec)]])
         )
    B = T @ A

    velocity_M = B @ np.array([
        [R_v],
        [k*mu_ra_cosdec*R],
        [k*mu_dec*R]
    ])

    return velocity_M

def get_rectangular_heliocentric_coords(l, b, R):
    return R * np.array([cos(b) * cos(l), cos(b)*sin(l), sin(b)])

def main():
    solar_velocity = [11.1, 247.2983, 7.25]
    solar_rect_coords = [-8.178, 0, 0.0208]
    data = pd.read_csv("orbits_table.txt", sep="\s+", on_bad_lines='warn')[["Cluster", "RA", "DEC", "Rsun", "<RV>", "mualpha", "mu_delta"]]
    choosed_clusters = pd.read_csv("choosed_clusters.csv")
    data = pd.merge(choosed_clusters, data, how="left", on="Cluster")
    data[["RA", "DEC", "Rsun", "<RV>", "mualpha", "mu_delta"]] = data[["RA", "DEC", "Rsun", "<RV>", "mualpha", "mu_delta"]].astype(float)
    T = get_T()
    for i in range(len(data)):
        data.loc[i, "l"], data.loc[i, "b"] = transform_coords(data.loc[i, "RA"], data.loc[i, "DEC"], T)
        data.loc[i, "U"], data.loc[i, "V"], data.loc[i, "W"] = get_velocity(
            data.loc[i, "RA"],
            data.loc[i, "DEC"],
            data.loc[i, "<RV>"],
            data.loc[i, "mualpha"],
            data.loc[i, "mu_delta"],
            data.loc[i, "Rsun"],
            T
        )
        data.loc[i, "x"], data.loc[i, "y"], data.loc[i, "z"] = get_rectangular_heliocentric_coords(
            data.loc[i, "l"],
            data.loc[i, "b"],
            data.loc[i, "Rsun"]
        )
    data[["U"]], data[["V"]], data[['W']] = \
        data[["U"]] + solar_velocity[0], data[["V"]] + solar_velocity[1], data[['W']] + solar_velocity[2]
    data[["x"]], data[["y"]], data[['z']] = \
        data[["x"]] + solar_rect_coords[0], data[["y"]] + solar_rect_coords[1], data[['z']] + solar_rect_coords[2]

    data.to_csv("my_orbits_table.csv")

    
main()