"""Vectorised element matrices.

T3: linear triangle (2D section).
W6: 6-node linear wedge (triangular prism extruded in z). As the wedges
come from a vertical extrusion, the Jacobian is separable: triangle in (x, y) times
segment in z. Integration: 3 points on the triangle x 2 Gauss points in z.
"""
import numpy as np

# Integration points: area coordinates on the triangle and zeta in z
GP_TRI = np.array([[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]])
GP_Z = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])


# ----------------------------------------------------------------------------
# Linear triangle (2D)
# ----------------------------------------------------------------------------
def geometry_t3(xy, tri):
    p = xy[tri]
    x, y = p[..., 0], p[..., 1]
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
    A2 = (x * b).sum(1)
    return 0.5 * A2, b / A2[:, None], c / A2[:, None]   # area, dN/dx, dN/dy


def conduction_t3(xy, tri):
    """Conduction matrix with k = 1 (per unit thickness) and area."""
    A, dx, dy = geometry_t3(xy, tri)
    K = A[:, None, None] * (dx[:, :, None] * dx[:, None, :] + dy[:, :, None] * dy[:, None, :])
    return K, A


# ----------------------------------------------------------------------------
# Linear wedge (3D)
# ----------------------------------------------------------------------------
def geometry_w6(xyz, conn3):
    A, dLdx, dLdy = geometry_t3(xyz[:, :2], conn3[:, :3])
    dz = xyz[conn3[:, 3], 2] - xyz[conn3[:, 0], 2]
    return A, dz, dLdx, dLdy


def _gradients_w6(A, dz, dLdx, dLdy, L, zeta):
    """Gradients of the 6 shape functions at a point (coordinates L, zeta):
    returns G (ne, 6, 3) and N (6,)."""
    fb, ft = (1 - zeta) / 2, (1 + zeta) / 2
    ne = len(A)
    G = np.empty((ne, 6, 3))
    G[:, :3, 0] = dLdx * fb; G[:, 3:, 0] = dLdx * ft
    G[:, :3, 1] = dLdy * fb; G[:, 3:, 1] = dLdy * ft
    G[:, :3, 2] = -L[None, :] / dz[:, None]
    G[:, 3:, 2] = L[None, :] / dz[:, None]
    N = np.concatenate([L * fb, L * ft])
    return G, N


def points_w6():
    for L in GP_TRI:
        for z in GP_Z:
            yield L, z


def conduction_w6(xyz, conn3):
    """Conduction matrix with k = 1 and volume of each wedge."""
    A, dz, dLdx, dLdy = geometry_w6(xyz, conn3)
    w = A * dz / 6.0                      # weight of each of the 6 points
    K = np.zeros((len(A), 6, 6))
    for L, z in points_w6():
        G, _ = _gradients_w6(A, dz, dLdx, dLdy, L, z)
        K += w[:, None, None] * np.einsum("eik,ejk->eij", G, G)
    return K, A * dz


def matrix_D(E, nu):
    """Isotropic 3D constitutive matrix, Voigt (xx, yy, zz, xy, yz, zx), per element."""
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
    """Strain-displacement matrix (ne, 6, 18) from the gradients."""
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


def stiffness_w6(xyz, conn3, E, nu, eps_th_nodes=None):
    """Stiffness (ne, 18, 18) and, if the nodal thermal strain per element is given
    (ne, 6), the vector of equivalent thermal forces (ne, 18)."""
    A, dz, dLdx, dLdy = geometry_w6(xyz, conn3)
    w = A * dz / 6.0
    D = matrix_D(E, nu)
    K = np.zeros((len(A), 18, 18))
    f = np.zeros((len(A), 18)) if eps_th_nodes is not None else None
    one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for L, z in points_w6():
        G, N = _gradients_w6(A, dz, dLdx, dLdy, L, z)
        B = _B_w6(G)
        DB = np.einsum("eij,ejk->eik", D, B)
        K += w[:, None, None] * np.einsum("eji,ejk->eik", B, DB)
        if f is not None:
            e_gp = eps_th_nodes @ N                          # (ne,)
            sig0 = (D @ one) * e_gp[:, None]                 # (ne, 6)
            f += w[:, None] * np.einsum("eji,ej->ei", B, sig0)
    return K, f


def stiffness_w6_D(xyz, conn3, D, eps_th_nodes=None):
    """As stiffness_w6, but with the constitutive matrix D (ne, 6, 6) given per
    element (orthotropic materials, such as the joints)."""
    A, dz, dLdx, dLdy = geometry_w6(xyz, conn3)
    w = A * dz / 6.0
    K = np.zeros((len(A), 18, 18))
    f = np.zeros((len(A), 18)) if eps_th_nodes is not None else None
    one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for L, z in points_w6():
        G, N = _gradients_w6(A, dz, dLdx, dLdy, L, z)
        B = _B_w6(G)
        DB = np.einsum("eij,ejk->eik", D, B)
        K += w[:, None, None] * np.einsum("eji,ejk->eik", B, DB)
        if f is not None:
            e_gp = eps_th_nodes @ N
            sig0 = (D @ one) * e_gp[:, None]
            f += w[:, None] * np.einsum("eji,ej->ei", B, sig0)
    return K, f


def thermal_force_w6(xyz, conn3, E, nu, eps_th_nodes):
    """Nodal forces equivalent to the thermal strain (ne, 18), isotropic material."""
    A, dz, dLdx, dLdy = geometry_w6(xyz, conn3)
    w = A * dz / 6.0
    D = matrix_D(E, nu)
    f = np.zeros((len(A), 18))
    one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for L, z in points_w6():
        G, N = _gradients_w6(A, dz, dLdx, dLdy, L, z)
        B = _B_w6(G)
        sig0 = (D @ one) * (eps_th_nodes @ N)[:, None]
        f += w[:, None] * np.einsum("eji,ej->ei", B, sig0)
    return f


def B_centroid_w6(xyz, conn3):
    """B matrix at the centroid (ne, 6, 18) and the N vector at the centroid (6,)."""
    A, dz, dLdx, dLdy = geometry_w6(xyz, conn3)
    G, N = _gradients_w6(A, dz, dLdx, dLdy, np.full(3, 1 / 3), 0.0)
    return _B_w6(G), N


def stresses_w6(xyz, conn3, E, nu, u_el, eps_th_nodes):
    """Stresses at the centroid (ne, 6) from the displacements (ne, 18)."""
    A, dz, dLdx, dLdy = geometry_w6(xyz, conn3)
    G, N = _gradients_w6(A, dz, dLdx, dLdy, np.full(3, 1 / 3), 0.0)
    B = _B_w6(G)
    eps = np.einsum("eij,ej->ei", B, u_el)
    e_c = eps_th_nodes @ N
    eps[:, :3] -= e_c[:, None]
    return np.einsum("eij,ej->ei", matrix_D(E, nu), eps)


def principal_stresses(sig):
    """Principal stresses (ne, 3) in increasing order, from Voigt (ne, 6)."""
    S = np.empty(sig.shape[:-1] + (3, 3))
    S[..., 0, 0], S[..., 1, 1], S[..., 2, 2] = sig[..., 0], sig[..., 1], sig[..., 2]
    S[..., 0, 1] = S[..., 1, 0] = sig[..., 3]
    S[..., 1, 2] = S[..., 2, 1] = sig[..., 4]
    S[..., 0, 2] = S[..., 2, 0] = sig[..., 5]
    return np.linalg.eigvalsh(S)
