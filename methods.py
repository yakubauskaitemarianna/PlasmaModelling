def Runge_Kutta_IV_3var_1step(f, phi, x0, y0, z0, h):
    k1 = f(x0, y0, z0)
    p1 = phi(x0, y0, z0)

    k2 = f(x0 + h/2, y0 + k1 * h / 2, z0 + p1 * h / 2)
    p2 = phi(x0 + h/2, y0 + k1 * h / 2, z0 + p1 * h / 2)

    k3 = f(x0 + h / 2, y0 + k2 * h / 2, z0 + p2 * h / 2)
    p3 = phi(x0 + h / 2, y0 + k2 * h / 2, z0 + p2 * h / 2)

    k4 = f(x0 + h / 2, y0 + k3 * h, z0 + p3 * h)
    p4 = phi(x0 + h / 2, y0 + k3 * h, z0 + p3 * h)

    zn = z0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    yn = y0 + h / 6 * (p1 + 2 * p2 + 2 * p3 + p4)

    return yn, zn



