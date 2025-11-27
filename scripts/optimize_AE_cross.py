import pyfurnace as pf
from scipy.optimize import minimize
import numpy as np


def crossover_symmetry_error(params, helix_length):
    ae_v_shift, ae_rot, ae_hx_shift, ae_hy_shift = params

    pf.Coords.set_AE_crossover_params(
        ae_v_shift=ae_v_shift,
        ae_rot=ae_rot * np.pi / 180,  # convert to radians
        ae_hx_shift=ae_hx_shift,
        ae_hy_shift=ae_hy_shift,
    )

    dae_T, dae_T_inv, dae_T_2, dae_T_inv_2 = pf.Coords.compute_AE_crossover(
        helix_length, use_cached=False, compute_all_T=True
    )

    # use Frobenius norm of matrix differences
    diff1 = dae_T - dae_T_2
    diff2 = dae_T_inv - dae_T_inv_2

    # you can optionally weight them, or add regularization terms if you want
    err = np.linalg.norm(diff1) ** 2 + np.linalg.norm(diff2) ** 2
    return err


helix_length = 11  # or whatever you use

# initial guess: your current values
x0 = np.array([0.99, 40.52, 2.23, -0.52])

# (optional) bounds if you want to keep parameters in a sane range
bounds = [
    (-2.0, 2.0),  # dae_v_shift
    (-90.0, 90.0),  # dae_rot
    (-5.0, 5.0),  # dae_hx_shift (radius-ish)
    (-5.0, 5.0),  # dae_hy_shift
]

res = minimize(
    crossover_symmetry_error,
    x0,
    args=(helix_length,),
    method="Nelder-Mead",  # or "L-BFGS-B" if you use 'bounds'
    # bounds=bounds,
    options={"maxiter": 200, "xatol": 1e-6, "fatol": 1e-8},
)

print("Optimization success:", res.success)
print("Message:", res.message)
print()
print("Optimized parameters:")
dae_v_opt, dae_rot_opt, dae_hx_opt, dae_hy_opt = res.x
print("  ae_v_shift =", dae_v_opt)
print("  ae_rot     =", dae_rot_opt)
print("  ae_hx_shift=", dae_hx_opt)
print("  ae_hy_shift=", dae_hy_opt)

print("Final symmetry error =", crossover_symmetry_error(res.x, helix_length))

print()
print("Final AE T matrices:")
pf.Coords.set_AE_crossover_params(
    ae_v_shift=dae_v_opt,
    ae_rot=dae_rot_opt,
    ae_hx_shift=dae_hx_opt,
    ae_hy_shift=dae_hy_opt,
)
dae_T, dae_T_inv = pf.Coords.compute_AE_crossover(
    helix_length, use_cached=False, compute_all_T=False
)
print("DAE T:\n", dae_T)
print("DAE T inv:\n", dae_T_inv)
