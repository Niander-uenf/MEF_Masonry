"""Matrizes de elemento vetorizadas.

T3: triângulo linear (seção 2D).
W6: cunha linear de 6 nós (prisma triangular extrudado em z). Como as cunhas
vêm de extrusão vertical, o jacobiano é separável: triângulo em (x, y) vezes
segmento em z. Integração: 3 pontos no triângulo x 2 pontos de Gauss em z.
"""
import numpy as np

# Pontos de integração: coordenadas de área no triângulo e zeta em z
GP_TRI = np.array([[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]])
GP_Z = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])


# ----------------------------------------------------------------------------
# Triângulo linear (2D)
# ----------------------------------------------------------------------------
def geometria_t3(xy, tri):
    p = xy[tri]
    x, y = p[..., 0], p[..., 1]
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
    A2 = (x * b).sum(1)
    return 0.5 * A2, b / A2[:, None], c / A2[:, None]   # área, dN/dx, dN/dy


def conducao_t3(xy, tri):
    """Matriz de condução com k = 1 (por unidade de espessura) e área."""
    A, dx, dy = geometria_t3(xy, tri)
    K = A[:, None, None] * (dx[:, :, None] * dx[:, None, :] + dy[:, :, None] * dy[:, None, :])
    return K, A


# ----------------------------------------------------------------------------
# Cunha linear (3D)
# ----------------------------------------------------------------------------
def geometria_w6(xyz, cun):
    A, dLdx, dLdy = geometria_t3(xyz[:, :2], cun[:, :3])
    dz = xyz[cun[:, 3], 2] - xyz[cun[:, 0], 2]
    return A, dz, dLdx, dLdy


def _gradientes_w6(A, dz, dLdx, dLdy, L, zeta):
    """Gradientes das 6 funções de forma num ponto (coordenadas L, zeta):
    devolve G (ne, 6, 3) e N (6,)."""
    fb, ft = (1 - zeta) / 2, (1 + zeta) / 2
    ne = len(A)
    G = np.empty((ne, 6, 3))
    G[:, :3, 0] = dLdx * fb; G[:, 3:, 0] = dLdx * ft
    G[:, :3, 1] = dLdy * fb; G[:, 3:, 1] = dLdy * ft
    G[:, :3, 2] = -L[None, :] / dz[:, None]
    G[:, 3:, 2] = L[None, :] / dz[:, None]
    N = np.concatenate([L * fb, L * ft])
    return G, N


def pontos_w6():
    for L in GP_TRI:
        for z in GP_Z:
            yield L, z


def conducao_w6(xyz, cun):
    """Matriz de condução com k = 1 e volume de cada cunha."""
    A, dz, dLdx, dLdy = geometria_w6(xyz, cun)
    w = A * dz / 6.0                      # peso de cada um dos 6 pontos
    K = np.zeros((len(A), 6, 6))
    for L, z in pontos_w6():
        G, _ = _gradientes_w6(A, dz, dLdx, dLdy, L, z)
        K += w[:, None, None] * np.einsum("eik,ejk->eij", G, G)
    return K, A * dz


def matriz_D(E, nu):
    """Matriz constitutiva isotrópica 3D, Voigt (xx, yy, zz, xy, yz, zx), por elemento."""
    E = np.asarray(E, float); nu = np.asarray(nu, float)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    D = np.zeros(E.shape + (6, 6))
    for i in range(3):
        for j in range(3):
            D[..., i, j] = lam
        D[..., i, i] = lam + 2 * mu
        D[..., 3 + i, 3 + i] = mu
    return D


def _B_w6(G):
    """Matriz deformação-deslocamento (ne, 6, 18) a partir dos gradientes."""
    ne = G.shape[0]
    B = np.zeros((ne, 6, 18))
    gx, gy, gz = G[..., 0], G[..., 1], G[..., 2]
    B[:, 0, 0::3] = gx
    B[:, 1, 1::3] = gy
    B[:, 2, 2::3] = gz
    B[:, 3, 0::3] = gy; B[:, 3, 1::3] = gx
    B[:, 4, 1::3] = gz; B[:, 4, 2::3] = gy
    B[:, 5, 0::3] = gz; B[:, 5, 2::3] = gx
    return B


def rigidez_w6(xyz, cun, E, nu, eps_th_nos=None):
    """Rigidez (ne, 18, 18) e, se dada a deformação térmica nodal por elemento
    (ne, 6), o vetor de forças térmicas equivalentes (ne, 18)."""
    A, dz, dLdx, dLdy = geometria_w6(xyz, cun)
    w = A * dz / 6.0
    D = matriz_D(E, nu)
    K = np.zeros((len(A), 18, 18))
    f = np.zeros((len(A), 18)) if eps_th_nos is not None else None
    um = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for L, z in pontos_w6():
        G, N = _gradientes_w6(A, dz, dLdx, dLdy, L, z)
        B = _B_w6(G)
        DB = np.einsum("eij,ejk->eik", D, B)
        K += w[:, None, None] * np.einsum("eji,ejk->eik", B, DB)
        if f is not None:
            e_gp = eps_th_nos @ N                          # (ne,)
            sig0 = (D @ um) * e_gp[:, None]                 # (ne, 6)
            f += w[:, None] * np.einsum("eji,ej->ei", B, sig0)
    return K, f


def rigidez_w6_D(xyz, cun, D, eps_th_nos=None):
    """Como rigidez_w6, mas com a matriz constitutiva D (ne, 6, 6) dada por
    elemento (materiais ortotrópicos, como as juntas)."""
    A, dz, dLdx, dLdy = geometria_w6(xyz, cun)
    w = A * dz / 6.0
    K = np.zeros((len(A), 18, 18))
    f = np.zeros((len(A), 18)) if eps_th_nos is not None else None
    um = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for L, z in pontos_w6():
        G, N = _gradientes_w6(A, dz, dLdx, dLdy, L, z)
        B = _B_w6(G)
        DB = np.einsum("eij,ejk->eik", D, B)
        K += w[:, None, None] * np.einsum("eji,ejk->eik", B, DB)
        if f is not None:
            e_gp = eps_th_nos @ N
            sig0 = (D @ um) * e_gp[:, None]
            f += w[:, None] * np.einsum("eji,ej->ei", B, sig0)
    return K, f


def forca_termica_w6(xyz, cun, E, nu, eps_th_nos):
    """Forças nodais equivalentes à deformação térmica (ne, 18), material isotrópico."""
    A, dz, dLdx, dLdy = geometria_w6(xyz, cun)
    w = A * dz / 6.0
    D = matriz_D(E, nu)
    f = np.zeros((len(A), 18))
    um = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for L, z in pontos_w6():
        G, N = _gradientes_w6(A, dz, dLdx, dLdy, L, z)
        B = _B_w6(G)
        sig0 = (D @ um) * (eps_th_nos @ N)[:, None]
        f += w[:, None] * np.einsum("eji,ej->ei", B, sig0)
    return f


def B_centroide_w6(xyz, cun):
    """Matriz B no centróide (ne, 6, 18) e o vetor N no centróide (6,)."""
    A, dz, dLdx, dLdy = geometria_w6(xyz, cun)
    G, N = _gradientes_w6(A, dz, dLdx, dLdy, np.full(3, 1 / 3), 0.0)
    return _B_w6(G), N


def tensoes_w6(xyz, cun, E, nu, u_el, eps_th_nos):
    """Tensões no centróide (ne, 6) a partir dos deslocamentos (ne, 18)."""
    A, dz, dLdx, dLdy = geometria_w6(xyz, cun)
    G, N = _gradientes_w6(A, dz, dLdx, dLdy, np.full(3, 1 / 3), 0.0)
    B = _B_w6(G)
    eps = np.einsum("eij,ej->ei", B, u_el)
    e_c = eps_th_nos @ N
    eps[:, :3] -= e_c[:, None]
    return np.einsum("eij,ej->ei", matriz_D(E, nu), eps)


def principais(sig):
    """Tensões principais (ne, 3) em ordem crescente, a partir de Voigt (ne, 6)."""
    S = np.empty(sig.shape[:-1] + (3, 3))
    S[..., 0, 0], S[..., 1, 1], S[..., 2, 2] = sig[..., 0], sig[..., 1], sig[..., 2]
    S[..., 0, 1] = S[..., 1, 0] = sig[..., 3]
    S[..., 1, 2] = S[..., 2, 1] = sig[..., 4]
    S[..., 0, 2] = S[..., 2, 0] = sig[..., 5]
    return np.linalg.eigvalsh(S)
