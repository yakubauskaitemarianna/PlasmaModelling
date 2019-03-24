class Data:
    t_0 = 0
    t_n = 2e-4
    l_e = 12        # Distance between lamp electrodes
    Rad = 0.35      # Radius of the tube (upper bound of the integral)
    p_0 = 0.5
    T_0 = 300
    T_w = 4000
    L_k = 60e-6     # Inductance
    C_k = 150e-6    # Capacitance capacity
    R_k = 0.5       # Resistance
    U_0 = 1500      # The capacitor voltage at the initial time
    I_0 = 0.5         # The current strength of the circuit at the initial time
    tau = 1e-6
    eps = 1e-4
    p = 15

    @staticmethod
    def get_table_I_T0():
        return [[0.5, 6400, 0.4],
                [1.0, 6790, 0.55],
                [5.0, 7150, 1.7],
                [10.0, 7270, 3],
                [50.0, 8010, 11],
                [200.0, 9185, 32],
                [400.0, 10010, 40],
                [800.0, 11140, 41],
                [1200.0, 12010, 39]]

    @staticmethod
    def get_table_sigma():
        return [[4000,  0.031],
                [5000, 0.27],
                [6000, 2.05],
                [7000, 6.06],
                [8000, 12],
                [9000, 19.9],
                [10000, 29.6],
                [11000, 41.1],
                [12000, 54.1],
                [13000, 67.7],
                [14000, 81.5]]