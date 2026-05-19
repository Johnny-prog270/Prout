# main_poisson_1d.py
import argparse
import numpy as np
import matplotlib.pyplot as plt

from gmsh_utils import (
    gmsh_init, gmsh_finalize, build_1d_mesh,
    prepare_quadrature_and_basis, get_jacobians, end_dofs_from_nodes
)
from stiffness_sol import assemble_stiffness_and_rhs, assemble_stiffness_and_rhs_ex3
from dirichlet import solve_dirichlet
from errors_sol import compute_L2_H1_errors, compute_energy_gap

from plot_utils import plot_fe_solution_high_order, setup_interactive_figure

def computeErrors1d(elemType, elemTags, elemNodeTags, U, order, u_exact, grad_exact):
    xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
    jac, det, coords = get_jacobians(elemType, xi)
    errL2, errH1s, errH1 = compute_L2_H1_errors(elemType, elemTags, elemNodeTags, U, xi, w, N, gN, jac, det, coords, u_exact=u_exact, grad_exact=grad_exact)
    return errL2, errH1s, errH1

def computeEnergyGap1d(elemType, elemTags, elemNodeTags, U, order, kappa, f, u_exact, grad_exact):
    xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
    jac, det, coords = get_jacobians(elemType, xi)
    energyExact = (-1.0)* compute_energy_gap(elemType, elemTags, elemNodeTags, np.zeros_like(U), kappa, f, xi, w, N, gN, jac, det, coords, u_exact=u_exact, grad_exact=grad_exact)
    energyFem = compute_energy_gap(elemType, elemTags, elemNodeTags, U, kappa, f, xi, w, N, gN, jac, det, coords, u_exact=lambda x: 0.0, grad_exact=lambda x: np.zeros(3))
    energy_gap = compute_energy_gap(elemType, elemTags, elemNodeTags, U, kappa, f, xi, w, N, gN, jac, det, coords, u_exact=u_exact, grad_exact=grad_exact)
    return energyExact, energyFem, energy_gap

def main_ex1(cl1, cl2, L, order, a, b):

    gmsh_init("poisson_1d")
    _, elemType, _, nodeCoords, elemTags, elemNodeTags = build_1d_mesh(L, cl1, cl2, order)

    xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
    jac, det, coords = get_jacobians(elemType, xi)

    if np.isclose(a, 0.0):
        def u_exact(x): return -(x[0]**2/2 + b/6*x[0]**3 - (1.0/2 + b/6)*x[0])
    else:
        def u_exact(x): return -(b/(4*a)*x[0]**2 + (1.0/a - b/(2*a**2))*x[0] + 1.0/a*(-(4*a + b*a - 2*b)/(4*a * np.log(np.abs(1+a))))*np.log(np.abs(1+a*x[0])))
    def kappa(x): return 1.0 + a*x[0]
    def f(x): return 1.0 + b*x[0]
    if np.isclose(a, 0.0):
        def grad_exact(x): return np.array([-(x[0] + b/2*x[0]**2 - (1.0/2 + b/6)), 0.0, 0.0])
    else :
        def grad_exact(x): return np.array([-(b/(2*a)*x[0] + (1.0/a - b/(2*a**2)) + 1.0/a*(-(4*a + b*a - 2*b)/(4*a * np.log(np.abs(1+a))))*a/(1+a*x[0])), 0.0, 0.0])

    K_lil, F = assemble_stiffness_and_rhs(elemTags, elemNodeTags, jac, det, coords, w, N, gN, kappa, f)
    K = K_lil.tocsr()

    left, right = end_dofs_from_nodes(nodeCoords)
    dir_dofs = [left, right]
    dir_vals = [u_exact([0.0]), u_exact([L])]

    U = solve_dirichlet(K, F, dir_dofs, dir_vals)

    errL2, errH1s, errH1 =  computeErrors1d(elemType, elemTags, elemNodeTags, U, order+7, u_exact, grad_exact)

    print(f"Errors for a={a}, b={b}: L2={errL2:.6e}, H1s={errH1s:.6e}, H1={errH1:.6e}")

    fig, ax = plt.subplots()
    plot_fe_solution_high_order(elemType, elemNodeTags, nodeCoords, U, M=200, show_nodes=True, ax=ax, label="FE solution")
    gmsh_finalize()

    x_plot = np.linspace(0.0, L, 200)
    u_plot = u_exact([x_plot])
    ax.plot(x_plot, u_plot, "k--", label="Exact solution")
    ax.legend()
    ax.set_title(f"a = {a}, b = {b}")
    plt.show()

    return errL2, errH1s, errH1

def main_ex2a(cl1, cl2, L, order, a=0.0, b=0.0):

    gmsh_init("poisson_1d")
    _, elemType, _, nodeCoords, elemTags, elemNodeTags = build_1d_mesh(L, cl1, cl2, order)

    xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
    jac, det, coords = get_jacobians(elemType, xi)

    if np.isclose(a, 0.0):
        def u_exact(x): return -(x[0]**2/2 + b/6*x[0]**3 - (1.0/2 + b/6)*x[0])
    else:
        def u_exact(x): return -(b/(4*a)*x[0]**2 + (1.0/a - b/(2*a**2))*x[0] + 1.0/a*(-(4*a + b*a - 2*b)/(4*a * np.log(np.abs(1+a))))*np.log(np.abs(1+a*x[0])))
    def kappa(x): return 1.0 + a*x[0]
    def f(x): return 1.0 + b*x[0]
    if np.isclose(a, 0.0):
        def grad_exact(x): return np.array([-(x[0] + b/2*x[0]**2 - (1.0/2 + b/6)), 0.0, 0.0])
    else :
        def grad_exact(x): return np.array([-(b/(2*a)*x[0] + (1.0/a - b/(2*a**2)) + 1.0/a*(-(4*a + b*a - 2*b)/(4*a * np.log(np.abs(1+a))))*a/(1+a*x[0])), 0.0, 0.0])

    K_lil, F = assemble_stiffness_and_rhs(elemTags, elemNodeTags, jac, det, coords, w, N, gN, kappa, f)
    K = K_lil.tocsr()

    left, right = end_dofs_from_nodes(nodeCoords)
    dir_dofs = [left, right]
    dir_vals = [u_exact([0.0]), u_exact([L])]

    U = solve_dirichlet(K, F, dir_dofs, dir_vals)

    energyExact, energyFem, energyGap = computeEnergyGap1d(elemType, elemTags, elemNodeTags, U, order+7, kappa, f, u_exact, grad_exact)
    print(f"Energy of exact solution: {energyExact:.6e}"
          f"\nEnergy of FE solution: {energyFem:.6e}"
          f"\nEnergy gap: {energyGap:.6e}")

    fig, ax = plt.subplots()
    plot_fe_solution_high_order(
        elemType, elemNodeTags, nodeCoords, U,
        M=200, show_nodes=True, ax=ax, label="FE solution"
    )
    gmsh_finalize()
    x_plot = np.linspace(0.0, L, 200)
    u_plot = u_exact([x_plot])
    ax.plot(x_plot, u_plot, "k--", label="Exact solution")
    ax.legend()
    plt.show()

    return None

def main_ex2b(cl1, cl2, L, order, a=0.0, b=0.0):
    if np.isclose(a, 0.0):
        def u_exact(x): return -(x[0]**2/2 + b/6*x[0]**3 - (1.0/2 + b/6)*x[0])
    else:
        def u_exact(x): return -(b/(4*a)*x[0]**2 + (1.0/a - b/(2*a**2))*x[0] + 1.0/a*(-(4*a + b*a - 2*b)/(4*a * np.log(np.abs(1+a))))*np.log(np.abs(1+a*x[0])))
    def kappa(x): return 1.0 + a*x[0]
    def f(x): return 1.0 + b*x[0]
    if np.isclose(a, 0.0):
        def grad_exact(x): return np.array([-(x[0] + b/2*x[0]**2 - (1.0/2 + b/6)), 0.0, 0.0])
    else :
        def grad_exact(x): return np.array([-(b/(2*a)*x[0] + (1.0/a - b/(2*a**2)) + 1.0/a*(-(4*a + b*a - 2*b)/(4*a * np.log(np.abs(1+a))))*a/(1+a*x[0])), 0.0, 0.0])

    fig, ax = plt.subplots()
    cl1 = args.cl1
    cl2 = args.cl2
    EnergyGap = []
    hs = []
    for iter in range(args.n):
        gmsh_init("poisson_1d")
        _, elemType, _, nodeCoords, elemTags, elemNodeTags = build_1d_mesh(L, cl1, cl2, order)

        xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
        jac, det, coords = get_jacobians(elemType, xi)

        K_lil, F = assemble_stiffness_and_rhs(elemTags, elemNodeTags, jac, det, coords, w, N, gN, kappa, f)
        K = K_lil.tocsr()

        left, right = end_dofs_from_nodes(nodeCoords)
        dir_dofs = [left, right]
        dir_vals = [u_exact([0.0]), u_exact([L])]

        U = solve_dirichlet(K, F, dir_dofs, dir_vals)

        # optional exact gradient

        _, _,  energyGap = computeEnergyGap1d(elemType, elemTags, elemNodeTags, U, order+7, kappa, f, u_exact, grad_exact)
        hs.append(cl1)
        EnergyGap.append(energyGap)
        cl1 = cl1 / 1.35
        cl2 = cl2 / 1.35
    _, C, rmse_log = plot_convergence_loglog(hs, EnergyGap, order, ax, label = f"order {order}")
    ax.set_title("Error in Energy")
    plt.show()

    return

def main_ex3(cl1, cl2, order, r0, u0):

    gmsh_init("poisson_1d")
    _, elemType, _, nodeCoords, elemTags, elemNodeTags = build_1d_mesh(1.0-r0, cl1, cl2, order)

    xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
    jac, det, coords = get_jacobians(elemType, xi)

    def u_exact(x): return ((1.0 - u0*np.sqrt(r0))*np.sqrt(x[0]+r0) + (u0 - np.sqrt(r0))*np.sqrt(r0)/(np.sqrt(x[0]+r0)))/(1.0-r0)
    def kappa(x): return 2.0*(x[0] + r0)
    def f(x): return 0.0
    def grad_exact(x): return np.array([None, 0.0, 0.0])

    K_lil, F = assemble_stiffness_and_rhs_ex3(elemTags, elemNodeTags, jac, det, coords, w, N, gN, kappa, f)
    K = K_lil.tocsr()

    left, right = end_dofs_from_nodes(nodeCoords)
    dir_dofs = [left, right]
    dir_vals = [u0, 1.0]

    U = solve_dirichlet(K, F, dir_dofs, dir_vals)

    errL2, errH1s, errH1 =  computeErrors1d(elemType, elemTags, elemNodeTags, U, order+7, u_exact, grad_exact)

    fig, ax = plt.subplots()
    plot_fe_solution_high_order(elemType, elemNodeTags, nodeCoords, U, M=200, show_nodes=True, ax=ax, label="FE solution")
    gmsh_finalize()

    x_plot = np.linspace(0.0, 1.0-r0, 200)
    u_plot = u_exact([x_plot])
    ax.plot(x_plot, u_plot, "k--", label="Exact solution")
    ax.set_title(f"r0={r0}, u0={u0}")
    ax.set_xlabel("r-r0")
    ax.legend()
    plt.show()

    return

def main_conv(cl1, cl2, L, order):

    MM = 5.0
    def u_exact(x): return np.sin(MM*np.pi * x[0])
    def kappa(x): return 1.0
    def f(x): return (MM*MM*np.pi*np.pi)*np.sin(MM*np.pi*x[0])
    def grad_exact(x): return np.array([MM*np.pi*np.cos(MM*np.pi*x[0]), 0.0, 0.0])

    allH1 = []
    allL2 = []
    allhs = []
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for O in range (order) :
        cl1 = args.cl1
        cl2 = args.cl2
        H1 = []
        L2 = []
        hs = []
        for iter in range(args.n) :
            gmsh_init("poisson_1d")
            _, elemType, _, nodeCoords, elemTags, elemNodeTags = build_1d_mesh(L, cl1, cl2, order)

            xi, w, N, gN = prepare_quadrature_and_basis(elemType, order)
            jac, det, coords = get_jacobians(elemType, xi)

            K_lil, F = assemble_stiffness_and_rhs(
                elemTags, elemNodeTags, jac, det, coords, w, N, gN, kappa, f
            )
            K = K_lil.tocsr()

            left, right = end_dofs_from_nodes(nodeCoords)
            dir_dofs = [left, right]
            dir_vals = [u_exact([0.0]), u_exact([L])]

            U = solve_dirichlet(K, F, dir_dofs, dir_vals)

            # optional exact gradient

            errL2, errH1s, errH1 =  computeErrors1d (elemType, elemTags, elemNodeTags, U, order+7, u_exact, grad_exact)
            hs .append(cl1)
            L2.append(errL2)
            H1.append(errH1s)
            cl1 = cl1 / 1.35
            cl2 = cl2 / 1.35
            gmsh_finalize()
        allH1.append(H1)
        allL2.append(L2)
        allhs.append(hs)
        _, C, rmse_log = plot_convergence_loglog(hs, H1, O+1, ax1, label = f"order {O+1}")
        _, C, rmse_log = plot_convergence_loglog(hs, L2, O+2, ax2, label = f"order {O+1}")
    ax1.set_title("Error in H1 seminorm")
    ax2.set_title("Error in L2 norm")
    plt.show()

    return


def plot_convergence_loglog(h, err, p, ax=None, label="error", fit_label=None, show=True):
    """
    Plot log-log convergence curve and a best-aligned reference line with user-given slope p.

    We fit: err_fit(h) = C * h^p
    In log-space: log(err) ≈ p*log(h) + log(C)
    Best log(C) (least squares) is mean(log(err) - p*log(h)).

    Parameters
    ----------
    h : array-like
        Mesh sizes
    err : array-like
        Errors (positive)
    p : float
        Desired slope (e.g. 1, 2, k+1, ...)
    ax : matplotlib axis or None
    label : str
        Label for numerical curve
    fit_label : str or None
        Label for fitted reference line
    show : bool
        plt.show() if True

    Returns
    -------
    ax, C, rmse_log
    """
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)

    # sort by h (optional but nicer)
    idx = np.argsort(h)
    h = h[idx]
    err = err[idx]

    logh = np.log(h)
    loge = np.log(err)

    # best intercept in log-space for fixed slope p
    logC = np.mean(loge - p * logh)
    C = float(np.exp(logC))

    # fitted curve
    err_fit = C * h**p

    # quality metric (optional)
    rmse_log = float(np.sqrt(np.mean((loge - (p*logh + logC))**2)))

    if ax is None:
        fig, ax = plt.subplots()

    ax.loglog(h, err, "o", label=label)
    if fit_label is None:
        fit_label = rf"fit: $C h^{{{p}}}$ (C={C:.2e})"
    ax.loglog(h, err_fit, "--", label=fit_label)

    ax.set_xlabel("h")
    ax.set_ylabel("error")
    ax.grid(True, which="both")
    ax.legend()

    return ax, C, rmse_log

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Poisson 1D with Gmsh high-order FE")
    parser.add_argument("-order", type=int, default=1)
    parser.add_argument("-cl1", type=float, default=0.1)
    parser.add_argument("-cl2", type=float, default=0.1)
    parser.add_argument("-n", type=int, default=10)
    parser.add_argument("-L", type=float, default=1.0)
    parser.add_argument("-Ex", type=int, default=1)
    args = parser.parse_args()
    order = args.order
    L = args.L
    cl1 = args.cl1
    cl2 = args.cl2
    
    if args.Ex == 1 :
        errL2, errH1s, errH1 = main_ex1(cl1, cl2, L, order, 0, 0)
        errL2, errH1s, errH1 = main_ex1(cl1, cl2, L, order, 10, 0)
        errL2, errH1s, errH1 = main_ex1(cl1, cl2, L, order, 0, 10)
        errL2, errH1s, errH1 = main_ex1(cl1, cl2, L, order, 10, 10)
    elif args.Ex == 2 :
        main_ex2a(cl1, cl2, L, order, a=10.0, b=0.0)
        main_ex2b(cl1, cl2, L, order, a=10.0, b=0.0)
    elif args.Ex == 3 :
        main_ex3(cl1, cl2, order, r0=0.01, u0=0.5)
        main_ex3(cl1/3, cl2, order, r0=0.01, u0=0.5)
        main_ex3(cl1, cl2, order, r0=0.00, u0=0.0)
        main_ex3(cl1/3, cl2, order, r0=0.00, u0=0.0)
    elif args.Ex == 4 :
        errL2, errH1s, errH1 = main_conv(cl1, cl2, L, order)
    else:
        raise ValueError(f"Unknown exercise number: {args.Ex}")
