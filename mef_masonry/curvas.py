"""Curvas de temperatura do gás (forno ou incêndio padrão)."""
import numpy as np


def rampa_patamar(T_inicial=30.0, taxa_C_min=5.0, T_alvo=600.0, patamar_s=1800.0):
    """Forno elétrico: rampa linear até T_alvo e patamar. Devolve (f(t), t_final, t_fim_rampa)."""
    t_rampa = (T_alvo - T_inicial) / (taxa_C_min / 60.0)

    def f(t):
        return float(T_inicial + min(t, t_rampa) * taxa_C_min / 60.0)
    return f, t_rampa + patamar_s, t_rampa


def iso834(T_inicial=20.0):
    """Incêndio-padrão ISO 834 / NBR 5628: T = T0 + 345 log10(8 t + 1), t em min."""
    def f(t):
        return float(T_inicial + 345.0 * np.log10(8.0 * t / 60.0 + 1.0))
    return f


def constante(T):
    return lambda t: float(T)


def rampa_patamar_resfriamento(T_inicial=30.0, taxa_C_min=5.0, T_alvo=600.0, patamar_s=1800.0,
                               taxa_resfriamento_C_min=2.0, patamar_final_s=3600.0):
    """Ensaio residual: rampa, patamar, resfriamento linear até T_inicial e
    espera final. Devolve (f(t), t_final, (t_fim_rampa, t_fim_patamar, t_fim_resfriamento))."""
    t1 = (T_alvo - T_inicial) / (taxa_C_min / 60.0)
    t2 = t1 + patamar_s
    t3 = t2 + (T_alvo - T_inicial) / (taxa_resfriamento_C_min / 60.0)

    def f(t):
        if t <= t1:
            return float(T_inicial + t * taxa_C_min / 60.0)
        if t <= t2:
            return float(T_alvo)
        if t <= t3:
            return float(T_alvo - (t - t2) * taxa_resfriamento_C_min / 60.0)
        return float(T_inicial)
    return f, t3 + patamar_final_s, (t1, t2, t3)
