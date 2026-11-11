"""Gas temperature curves (furnace or standard fire)."""
import numpy as np


def ramp_dwell(T_initial=30.0, rate_C_min=5.0, T_target=600.0, dwell_s=1800.0):
    """Electric furnace: linear ramp to T_target and dwell. Returns (f(t), t_final, t_end_ramp)."""
    t_ramp = (T_target - T_initial) / (rate_C_min / 60.0)

    def f(t):
        return float(T_initial + min(t, t_ramp) * rate_C_min / 60.0)
    return f, t_ramp + dwell_s, t_ramp


def iso834(T_initial=20.0):
    """ISO 834 / NBR 5628 standard fire: T = T0 + 345 log10(8 t + 1), t in min."""
    def f(t):
        return float(T_initial + 345.0 * np.log10(8.0 * t / 60.0 + 1.0))
    return f


def constant(T):
    return lambda t: float(T)


def ramp_dwell_cooling(T_initial=30.0, rate_C_min=5.0, T_target=600.0, dwell_s=1800.0,
                               cooling_rate_C_min=2.0, final_dwell_s=3600.0):
    """Residual test: ramp, dwell, linear cooling to T_initial and a
    final wait. Returns (f(t), t_final, (t_end_ramp, t_end_dwell, t_end_cooling))."""
    t1 = (T_target - T_initial) / (rate_C_min / 60.0)
    t2 = t1 + dwell_s
    t3 = t2 + (T_target - T_initial) / (cooling_rate_C_min / 60.0)

    def f(t):
        if t <= t1:
            return float(T_initial + t * rate_C_min / 60.0)
        if t <= t2:
            return float(T_target)
        if t <= t3:
            return float(T_target - (t - t2) * cooling_rate_C_min / 60.0)
        return float(T_initial)
    return f, t3 + final_dwell_s, (t1, t2, t3)
