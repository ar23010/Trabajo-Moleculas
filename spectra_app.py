import math
from matplotlib.patches import Circle
import numpy as np
import plotly.graph_objects as go
import plotly.offline as pyo
import streamlit as st
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter

#----- Gráfico del Modelo 3D NH3 -----
def dibujar_nh3():
    # ---------- Parámetros ----------
    NH_LENGTH = 1.01
    SHADOW_Z  = -0.55
    COLORS = {
        'N': '#1B5E20',
        'H': '#A5D6A7',
        'BOND': '#66BB6A'
    }
    RADII = {'N': 0.34, 'H': 0.18, 'BOND': 0.10}

    # ---------- Utilidades geométricas ----------
    def sphere_mesh(center, r=1.0, nu=36, nv=18):
        u = np.linspace(0, 2*np.pi, nu)
        v = np.linspace(0, np.pi, nv)
        u, v = np.meshgrid(u, v)
        x = center[0] + r*np.cos(u)*np.sin(v)
        y = center[1] + r*np.sin(u)*np.sin(v)
        z = center[2] + r*np.cos(v)
        return x.ravel(), y.ravel(), z.ravel()

    def cylinder_mesh(p0, p1, r=0.10, seg=36, slc=12):
        p0 = np.asarray(p0, float); p1 = np.asarray(p1, float)
        v = p1 - p0
        L = np.linalg.norm(v)
        if L == 0:
            return np.array([]), np.array([]), np.array([])
        v = v / L

        a = np.array([1,0,0]) if abs(v[0]) < 0.9 else np.array([0,1,0])
        n1 = np.cross(v, a); n1 /= np.linalg.norm(n1)
        n2 = np.cross(v, n1); n2 /= np.linalg.norm(n2)

        th = np.linspace(0, 2*np.pi, seg)
        zz = np.linspace(0, L, slc)
        th, zz = np.meshgrid(th, zz)

        x = p0[0] + v[0]*zz + r*(n1[0]*np.cos(th) + n2[0]*np.sin(th))
        y = p0[1] + v[1]*zz + r*(n1[1]*np.cos(th) + n2[1]*np.sin(th))
        z = p0[2] + v[2]*zz + r*(n1[2]*np.cos(th) + n2[2]*np.sin(th))
        return x.ravel(), y.ravel(), z.ravel()

    # ---------- Piezas de la escena ----------
    def add_shadow(fig, c, r, darkness=0.35):
        t = np.linspace(0, 2*np.pi, 60)
        x = c[0] + (r*1.25)*np.cos(t)
        y = c[1] + (r*0.85)*np.sin(t)
        z = np.full_like(x, SHADOW_Z)
        fig.add_trace(go.Mesh3d(
            x=x, y=y, z=z,
            color=f'rgba(10,10,10,{darkness:.3f})',
            opacity=1.0, showscale=False, hoverinfo='skip',
            i=None, j=None, k=None
        ))

    def add_atom(fig, pos, kind):
        x, y, z = sphere_mesh(pos, RADII[kind], 40, 22)
        fig.add_trace(go.Mesh3d(
            x=x, y=y, z=z,
            color=COLORS[kind],
            opacity=1.0, alphahull=0, flatshading=False,
            name=kind, hovertext=kind, hoverinfo='text'
        ))

    def add_bond(fig, p0, p1):
        x, y, z = cylinder_mesh(p0, p1, RADII['BOND'])
        fig.add_trace(go.Mesh3d(
            x=x, y=y, z=z,
            color=COLORS['BOND'],
            opacity=1.0, name='bond', hoverinfo='skip'
        ))

    def add_plane(fig, size=2.6):
        xx = np.linspace(-size, size, 2)
        yy = np.linspace(-size, size, 2)
        XX, YY = np.meshgrid(xx, yy)
        ZZ = np.full_like(XX, SHADOW_Z - 0.002)
        fig.add_trace(go.Surface(
            x=XX, y=YY, z=ZZ,
            showscale=False, opacity=0.08,
            colorscale=[[0,'rgb(200,200,200)'],[1,'rgb(200,200,200)']],
            hoverinfo='skip'
        ))

    # ---------- Geometría NH3 ----------
    base_dirs = np.array([
        [ 1, -1, -1],
        [-1,  1, -1],
        [-1, -1,  1],
    ], dtype=float)
    base_dirs = base_dirs / np.linalg.norm(base_dirs, axis=1)[:, None] * NH_LENGTH

    N = np.array([0.0, 0.0, 0.0])
    H1, H2, H3 = base_dirs

    # ---------- Construcción de la figura ----------
    fig = go.Figure()

    add_shadow(fig, N,  RADII['N'], darkness=0.45)
    for H in (H1, H2, H3):
        add_shadow(fig, H, RADII['H'], darkness=0.25)

    add_atom(fig, N,  'N')
    add_atom(fig, H1, 'H')
    add_atom(fig, H2, 'H')
    add_atom(fig, H3, 'H')

    add_bond(fig, N, H1)
    add_bond(fig, N, H2)
    add_bond(fig, N, H3)

    add_plane(fig)

    fig.update_layout(
        title="Amoniaco (NH₃) — vista 3D (tonos verdes)",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode='data',
            camera=dict(eye=dict(x=1.6, y=-1.2, z=0.9))
        ),
        margin=dict(l=10, r=10, t=50, b=10),
        showlegend=False,
    )
    return fig

#----- Espectro IR NH3 -----
def dibujar_ir():
    # ---------- Eje de número de onda (4000 → 400 cm⁻1) ----------
    wavenumber = np.linspace(4000, 400, 2000)

    # ---------- Función generadora de picos gaussianos ----------
    def gauss(x, centro, ancho, intensidad):
        return intensidad * np.exp(-0.5 * ((x - centro) / ancho) ** 2)

    # ---------- Señales simuladas ----------
    polyamide6 = (
        gauss(wavenumber, 3300, 80, 1.0) +
        gauss(wavenumber, 1630, 50, 0.8) +
        gauss(wavenumber, 1540, 40, 0.6)
    )

    water = (
        gauss(wavenumber, 3400, 120, 1.0) +
        gauss(wavenumber, 1650, 80, 0.7)
    )

    mix = polyamide6 + 0.8 * water

    # ---------- Gráfico ----------
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(wavenumber, polyamide6, label="Poliamida 6", color="black", lw=1.4)
    ax.plot(wavenumber, water, label="Agua", color="green", lw=1.4)
    ax.plot(wavenumber, mix, label="Poliamida 6 + Agua", color="red", lw=1.4)

    ax.invert_xaxis()

    ax.set_title("Espectro IR simulado")
    ax.set_xlabel("Número de onda (cm⁻¹)")
    ax.set_ylabel("Absorbancia (a.u.)")
    ax.legend()
    ax.grid(alpha=0.3)

    return fig

#----- Gráfico Optimizado -----
def dibujar_conjunto_nh3():

    # --- Apariencia  ---
    COL_N = "#1B5E20"
    COL_H = "#A5D6A7"
    COL_BOND = "#81C784"

    # --- N (centros) a lo largo de una cadena zig-zag---
    main_atoms = np.array([
        [0.0, 6.0, 0.0],
        [0.3, 5.0, 0.0],
        [-0.3, 4.0, 0.0],
        [0.3, 3.0, 0.0],
        [-0.3, 2.0, 0.0],
        [0.3, 1.0, 0.0],
        [-0.3, 0.0, 0.0]
    ], dtype=float)

    # ---------- Helpers geométricos ----------
    def create_sphere(center, radius, color):
        u, v = np.mgrid[0:2*np.pi:12j, 0:np.pi:6j]
        x = radius * np.cos(u) * np.sin(v) + center[0]
        y = radius * np.sin(u) * np.sin(v) + center[1]
        z = radius * np.cos(v) + center[2]
        return go.Surface(
            x=x, y=y, z=z,
            colorscale=[[0, color], [1, color]],
            showscale=False, hoverinfo='skip', opacity=1
        )

    def create_cylinder(p0, p1, radius, color, resolution=8):
        p0, p1 = np.array(p0, float), np.array(p1, float)
        v = p1 - p0
        L = np.linalg.norm(v)
        if L == 0:
            return None
        v /= L

        a = np.array([1, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1, 0])
        n1 = np.cross(v, a); n1 /= np.linalg.norm(n1)
        n2 = np.cross(v, n1); n2 /= np.linalg.norm(n2)

        t = np.linspace(0, L, 2)
        th = np.linspace(0, 2*np.pi, resolution)
        t, th = np.meshgrid(t, th)

        X = p0[0] + v[0]*t + radius*np.sin(th)*n1[0] + radius*np.cos(th)*n2[0]
        Y = p0[1] + v[1]*t + radius*np.sin(th)*n1[1] + radius*np.cos(th)*n2[1]
        Z = p0[2] + v[2]*t + radius*np.sin(th)*n1[2] + radius*np.cos(th)*n2[2]
        return go.Surface(
            x=X, y=Y, z=Z,
            colorscale=[[0, color], [1, color]],
            showscale=False, hoverinfo='skip'
        )

    # ---------- Offsets de H para NH3 ----------
    def nh3_h_offsets(bond_len=0.55, lift=0.15):
        base = np.array([
            [ 1,  1,  1],
            [ 1, -1, -1],
            [-1,  1, -1],
        ], dtype=float)
        base = base / np.linalg.norm(base, axis=1)[:, None]
        vecs = base * bond_len
        vecs[:, 1] += lift
        return vecs.tolist()

    # ---------- Construcción de la figura ----------
    fig = go.Figure()

    # N (esferas) y enlaces N–N entre consecutivos
    for i, pos in enumerate(main_atoms):
        fig.add_trace(create_sphere(pos, radius=0.25, color=COL_N))
        if i < len(main_atoms) - 1:
            fig.add_trace(create_cylinder(main_atoms[i], main_atoms[i+1], 0.04, COL_BOND))

    # Para cada N, añadimos 3 H (satélites) con enlaces N–H
    for center in main_atoms:
        for off in nh3_h_offsets():
            h_pos = center + np.array(off, float)
            fig.add_trace(create_sphere(h_pos, radius=0.12, color=COL_H))
            fig.add_trace(create_cylinder(center, h_pos, 0.025, COL_BOND))

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode='data',
            bgcolor='rgb(14,17,23)'
        ),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    return fig

#----- Gráfico del Modelo 2D NH3 -----
def dibujar_nh3_2d():
    # ---- Apariencia ----
    COL_N_FILL = "#1B5E20"
    COL_N_EDGE = "#0D3D0D"

    COL_H_FILL = "#A5D6A7"
    COL_H_EDGE = "#4CAF50"

    COL_BOND = "#66BB6A"
    COL_BOND_INNER = "#C8E6C9"
    BG_COLOR = "#FFFFFF"

    R_N, R_H = 1.45, 1.05
    BOND_W, BOND_W_INNER = 28, 12

    # ---- Geometría NH3 ----
    NH = 3.0
    ANGLE = 120
    N = np.array([0.0, 0.0])

    H_positions = []
    for i in range(3):
        theta = math.radians(90 + i * ANGLE)
        H_positions.append(N + np.array([NH * math.cos(theta),
                                         -NH * math.sin(theta)]))

    # ---- Helpers ----
    def draw_bond(ax, p0, p1):
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                color=COL_BOND, linewidth=BOND_W, solid_capstyle="round", zorder=1)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                color=COL_BOND_INNER, linewidth=BOND_W_INNER, solid_capstyle="round", zorder=1.1)

    def glossy_sphere(ax, center, r, face, edge, edge_w=3.0,
                      highlight_shift=(-0.45, 0.45), highlight_alpha=0.40):
        shadow = Circle((center[0]+0.25, center[1]-0.25), r*1.02,
                        facecolor="black", edgecolor="none", alpha=0.18, zorder=0.4)
        ax.add_patch(shadow)
        base = Circle(center, r, facecolor=face, edgecolor=edge,
                      linewidth=edge_w, zorder=2)
        ax.add_patch(base)
        hx, hy = highlight_shift
        hi1 = Circle((center[0]+hx, center[1]+hy), r*0.78,
                     facecolor="white", edgecolor="none",
                     alpha=highlight_alpha*0.55, zorder=3)
        ax.add_patch(hi1)
        hi2 = Circle((center[0]+hx*0.55, center[1]+hy*0.55), r*0.38,
                     facecolor="white", edgecolor="none",
                     alpha=highlight_alpha, zorder=3)
        ax.add_patch(hi2)

    # ---- Render principal ----
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=160)
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    for H in H_positions:
        draw_bond(ax, N, H)

    glossy_sphere(ax, N,  R_N, COL_N_FILL, COL_N_EDGE,
                  edge_w=3.5, highlight_shift=(-0.55, 0.55))
    for H in H_positions:
        glossy_sphere(ax, H, R_H, COL_H_FILL, COL_H_EDGE,
                      highlight_shift=(-0.55, 0.55))

    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    ax.set_xlim(-6.5, 6.5)
    ax.set_ylim(-6.0, 6.0)

    fig.tight_layout()
    return fig

#----- Contenedor de Moléculas -----
def contenedor_molecular(
    n_molecules: int = 400,
    atoms_per_molecule: int = 6,
    box: float = 80.0,
    jitter_scale: float = 2.0,
    seed: int = 0,
    marker_size: int = 3,
    marker_opacity: float = 0.8,
    marker_color: str = "lightgreen",
):

    # ---------- Generación de moléculas ----------
    rng = np.random.default_rng(seed)

    def _random_molecule(center: np.ndarray) -> np.ndarray:
        """
        Crea 'atoms_per_molecule' puntos con distribución normal
        alrededor de 'center' con desviación 'jitter_scale'.
        """
        return rng.normal(0.0, jitter_scale, size=(atoms_per_molecule, 3)) + center

    centers = rng.uniform(2.0, box - 2.0, size=(n_molecules, 3))
    molecules = [_random_molecule(c) for c in centers]

    # ---------- Preparar datos para graficar ----------
    if molecules:
        coords = np.vstack(molecules)
        xs, ys, zs = coords[:, 0], coords[:, 1], coords[:, 2]
    else:
        xs = ys = zs = np.array([])

    # ---------- Construcción de la caja (aristas) ----------
    def _box_edges(L: float):
        edges = [
            (0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0), (0, 0, 0),

            (0, 0, L), (L, 0, L), (L, L, L), (0, L, L), (0, 0, L),

            (L, 0, L), (L, 0, 0), (L, L, 0), (L, L, L), (0, L, L), (0, L, 0)
        ]
        x, y, z = zip(*edges)
        return list(x), list(y), list(z)

    bx, by, bz = _box_edges(box)

    # ---------- Visualización 3D ----------
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs,
        mode="markers",
        marker=dict(size=marker_size, color=marker_color, opacity=marker_opacity),
        name="Átomos"
    ))

    fig.add_trace(go.Scatter3d(
        x=bx, y=by, z=bz,
        mode="lines",
        line=dict(width=3),
        name="Caja"
    ))

    fig.update_scenes(aspectmode="data")
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        showlegend=False,
        margin=dict(l=0, r=0, t=30, b=0),
        title=f"Conjunto de {n_molecules} moléculas sintéticas en caja {box} Å"
    )

    return fig

#----- Trabajo de Adhesión -----
def graficar_trabajo_adhesion():
    # ---------- Datos ----------
    x_h2o = np.array([0, 5, 10, 15, 20], dtype=float)

    y_vdw_contacto = np.array([2.1, 0.85, 0.60, 0.08, 0.05], dtype=float)
    y_ener_pot     = np.array([2.35, 1.05, 0.70, 0.10, 0.06], dtype=float)

    # ---------- Gráfico ----------
    fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=120)

    ax.plot(
        x_h2o, y_vdw_contacto,
        color="#FF02F2", marker="o", linewidth=2, label="VDW + Contacto"
    )
    ax.plot(
        x_h2o, y_ener_pot,
        color="#33FF00", marker="^", linewidth=2, label="Energía Potencial"
    )

    ax.set_title("Trabajo de Adhesión vs %NH₃", fontsize=13, fontweight="bold")
    ax.set_xlabel("% NH₃", fontsize=11)
    ax.set_ylabel("Energía (J/m²)", fontsize=11)

    ax.legend(title="Tipo de Energía")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    ax.set_xticks([0, 5, 10, 15, 20])

    fig.tight_layout()
    return fig

#----- Molécula Teórica -----
def mostrar_rdf():
    # ---------- Parámetros ----------
    n_bins = 75
    box_size = 3.0
    r_max = box_size / 2
    r_min = 0.25 

    # Carpeta de modelos
    ruta_modelos = os.path.join(os.path.dirname(__file__), "modelos")

    def leer_coordenadas_nh3(file_path: str) -> np.ndarray:
        coords = []
        try:
            with open(file_path, "r") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 6 and parts[0].endswith("NH3"):
                        x, y, z = map(float, parts[3:6])
                        coords.append([x, y, z])
        except Exception as e:
            st.warning(f"No se pudo leer {file_path}: {e}")
        return np.array(coords, dtype=float)

    def calcular_rdf(posiciones: np.ndarray, n_bins: int, r_max: float,
                     box_size: float, r_min: float = 0.25):
        
        if posiciones.size == 0:
            return np.array([]), np.array([])

        n = len(posiciones)
        distancias = []

        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(posiciones[i] - posiciones[j])
                if d >= r_min:
                    distancias.append(d)

        distancias = np.array(distancias, dtype=float)
        hist, edges = np.histogram(distancias, bins=n_bins, range=(0, r_max))

        r = 0.5 * (edges[1:] + edges[:-1])
        densidad = n / (box_size ** 3)

        shell_vol = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
        normalizacion = shell_vol * densidad

        normalizacion = np.where(normalizacion == 0, 1.0, normalizacion)
        g_r = hist / (normalizacion * n)
        return r, g_r


    n_puntos_demo = 200
    pos_random = np.random.rand(n_puntos_demo, 3) * box_size
    r_rand, g_rand = calcular_rdf(pos_random, n_bins, r_max, box_size, r_min=r_min)

    # Lista de modelos
    modelos = {
        "Modelo 1 (nh3.txt)": os.path.join(ruta_modelos, "nh3.txt"),
        "Modelo 2 (nh3_frame.txt)": os.path.join(ruta_modelos, "nh3_frame.txt"),
    }

    # ---------- Gráfico ----------
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=120)

    verde_base = "#F81313"
    verde_claro = "#CEB00A"
    verde_med   = "#016FFF"

    if r_rand.size:
        ax.plot(r_rand, g_rand, label="NH₃ (modelo teórico)", color=verde_claro, lw=2)

    ymax_vals = [g_rand.max() if g_rand.size else 1.0]

    for i, (nombre, path) in enumerate(modelos.items(), start=1):
        if os.path.exists(path):
            coords = leer_coordenadas_nh3(path)
            r, g = calcular_rdf(coords, n_bins, r_max, box_size, r_min=r_min)
            if r.size:
                color = [verde_base, verde_med][(i - 1) % 2]
                ax.plot(r, g, label=nombre, lw=2, color=color)
                ymax_vals.append(g.max() if g.size else 1.0)
            else:
                st.warning(f"Archivo sin datos válidos: {path}")
        else:
            st.warning(f"No se encontró el archivo: {path}")

    ax.set_xlabel("r (nm)")
    ax.set_ylabel("g(r)")
    ax.set_title("Función de Distribución Radial (RDF) — NH₃")
    ax.legend()
    ax.grid(alpha=0.35, linestyle="--", linewidth=0.8)

    ax.set_xlim(0, 1.5)
    ax.set_ylim(0, max(ymax_vals) * 1.10 if ymax_vals else 2.0)

    st.pyplot(fig)

#----- Raman -----
def mostrar_ir_raman(ruta_paso_1):
    """
    Lee un archivo de espectros (IR y Raman) y grafica ambos junto con tablas comparativas.
    """
    # --- Parsear archivo ---
    ir_data, raman_data = [], []
    seccion = None
    with open(ruta_paso_1, "r") as f:
        for linea in f:
            linea = linea.strip()
            if not linea or linea.startswith("#"):
                continue
            if linea.startswith("IR SPECTRUM"):
                seccion = "IR"
                continue
            elif linea.startswith("RAMAN SPECTRUM"):
                seccion = "RAMAN"
                continue
            elif linea[0].isdigit():  # Datos
                partes = linea.replace(":", "").replace("(", "").replace(")", "").split()
                if seccion == "IR":
                    modo, freq, eps, intensidad, t2, tx, ty, tz = partes
                    ir_data.append({
                        "Mode": int(modo),
                        "Freq (cm^-1)": float(freq),
                        "Eps (L/mol*cm)": float(eps),
                        "Int (km/mol)": float(intensidad),
                        "T^2 (a.u.)": float(t2),
                        "Tx": float(tx),
                        "Ty": float(ty),
                        "Tz": float(tz)
                    })
                elif seccion == "RAMAN":
                    modo, freq, actividad, depol = partes
                    raman_data.append({
                        "Mode": int(modo),
                        "Freq (cm^-1)": float(freq),
                        "Actividad": float(actividad),
                        "Depolarización": float(depol)
                    })

    # Convertir a DataFrames
    df_ir = pd.DataFrame(ir_data)
    df_raman = pd.DataFrame(raman_data)

    # --- Graficar IR ---
    fig_ir, ax_ir = plt.subplots(figsize=(7, 4))
    ax_ir.plot(df_ir["Freq (cm^-1)"], df_ir["Int (km/mol)"], "-o", color="blue")
    ax_ir.set_title("Espectro IR")
    ax_ir.set_xlabel("Frecuencia (cm⁻¹)")
    ax_ir.set_ylabel("Intensidad (km/mol)")

    # --- Graficar Raman ---
    fig_raman, ax_raman = plt.subplots(figsize=(7, 4))
    ax_raman.plot(df_raman["Freq (cm^-1)"], df_raman["Actividad"], "-o", color="green")
    ax_raman.set_title("Espectro Raman")
    ax_raman.set_xlabel("Frecuencia (cm⁻¹)")
    ax_raman.set_ylabel("Actividad")

    # --- Mostrar en Streamlit ---
    st.subheader("Tabla IR")
    st.dataframe(df_ir)

    st.subheader("Espectro IR")
    st.pyplot(fig_ir)

    st.subheader("Tabla Raman")
    st.dataframe(df_raman)

    st.subheader("Espectro Raman")
    st.pyplot(fig_raman)


# ----------------- Utilidades compartidas -----------------
_RE_XYZ = re.compile(
    r'^\s*([A-Z][a-z]?)\s+'
    r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+'
    r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+'
    r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$'
)
MASAS = {"H":1.00784, "C":12.0107, "N":14.0067, "O":15.999, "F":18.998, "P":30.9738, "S":32.065, "Cl":35.453, "Br":79.904, "I":126.904}

def _parse_orca_xyz(ruta: str):
    """Devuelve (elems:list[str], coords:np.ndarray[n,3])."""
    elems, coords, in_block = [], [], False
    with open(ruta, "r") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("*"):
                if not in_block and line.lower().startswith("* xyz"):
                    in_block = True
                    continue
                elif in_block:
                    break
                else:
                    continue
            if in_block:
                m = _RE_XYZ.match(line)
                if m:
                    elems.append(m.group(1))
                    coords.append([float(m.group(2)), float(m.group(3)), float(m.group(4))])
    return elems, (np.array(coords, float) if coords else np.empty((0,3)))

def _df_coordenadas(nombre: str, elems, coords) -> pd.DataFrame:
    if coords.size == 0:
        return pd.DataFrame(columns=["Molécula","Átomo","X(Å)","Y(Å)","Z(Å)"])
    df = pd.DataFrame({
        "Molécula": nombre,
        "Átomo": elems,
        "X(Å)": coords[:,0],
        "Y(Å)": coords[:,1],
        "Z(Å)": coords[:,2],
    })
    return df

def _fila_metricas(nombre: str, elems, coords) -> dict:
    n = len(elems)
    if n == 0:
        return {"Molécula": nombre, "Átomos": 0, "Masa (u)": 0.0, "Elementos": "-",
                "Centro (Å)": "-", "Caja Δx,Δy,Δz (Å)": "-", "R_g (Å)": 0.0,
                "Distancia media (Å)": 0.0, "Distancia máx (Å)": 0.0,
                "⟨N–H⟩ (Å)": "-", "⟨C–H⟩ (Å)": "-"}
    coords = np.asarray(coords, float)
    masa = float(sum(MASAS.get(el, 0.0) for el in elems))
    c = Counter(elems)
    comp = " ".join(f"{el}{c[el]}" for el in sorted(c))
    cen = coords.mean(axis=0); cen_s = f"({cen[0]:.3f}, {cen[1]:.3f}, {cen[2]:.3f})"
    dmin, dmax = coords.min(axis=0), coords.max(axis=0)
    box = dmax - dmin; box_s = f"({box[0]:.3f}, {box[1]:.3f}, {box[2]:.3f})"
    rg = float(np.sqrt(((coords - cen)**2).sum(axis=1).mean()))
    if n >= 2:
        i, j = np.triu_indices(n, k=1)
        dist = np.linalg.norm(coords[i]-coords[j], axis=1)
        dmed, dmx = float(dist.mean()), float(dist.max())
    else:
        dmed = dmx = 0.0

    # promedios geométricos de enlaces (umbral simple)
    def prom(el_a, el_b, rmin, rmax):
        idx_a = [k for k,e in enumerate(elems) if e==el_a]
        idx_b = [k for k,e in enumerate(elems) if e==el_b]
        vals=[]
        for ia in idx_a:
            for ib in idx_b:
                r = float(np.linalg.norm(coords[ia]-coords[ib]))
                if rmin <= r <= rmax:
                    vals.append(r)
        return f"{np.mean(vals):.3f}" if vals else "-"
    nh = prom("N","H",0.90,1.20)
    ch = prom("C","H",0.95,1.20)

    return {"Molécula": nombre, "Átomos": n, "Masa (u)": round(masa,4),
            "Elementos": comp, "Centro (Å)": cen_s, "Caja Δx,Δy,Δz (Å)": box_s,
            "R_g (Å)": round(rg,4), "Distancia media (Å)": round(dmed,4),
            "Distancia máx (Å)": round(dmx,4), "⟨N–H⟩ (Å)": nh, "⟨C–H⟩ (Å)": ch}

def _nh3_ref():
    """NH3 de referencia: N en (0,0,0), 3 H tetraédricos (~1.01 Å)."""
    NH = 1.01
    N = np.array([[0.0,0.0,0.0]])
    Hs = np.array([[ 1,-1,-1],[-1, 1,-1],[-1,-1, 1]], float)
    Hs = Hs/np.linalg.norm(Hs,axis=1)[:,None]*NH
    elems = ["N","H","H","H"]
    coords = np.vstack([N,Hs])
    return elems, coords

# ----- Comparación general ORCA vs NH3 -----
def comparar_moleculas_orca_vs_nh3(ruta: str) -> pd.DataFrame:
    st.header("Comparación · Molécula ORCA vs NH₃")

    if not os.path.exists(ruta):
        st.error(f"No se encontró el archivo:\n{ruta}")
        return pd.DataFrame()

    # Parsear
    elems_o, coords_o = _parse_orca_xyz(ruta)
    elems_ref, coords_ref = _nh3_ref()

    # Tablas de coordenadas
    st.subheader("Coordenadas ORCA")
    st.dataframe(_df_coordenadas("ORCA", elems_o, coords_o), use_container_width=True)

    st.subheader("Coordenadas NH₃ (referencia)")
    st.dataframe(_df_coordenadas("NH₃", elems_ref, coords_ref), use_container_width=True)

    # Gráfica
    fig, ax = plt.subplots(figsize=(7,4))
    if coords_o.size>0:
        ax.scatter(coords_o[:,0], coords_o[:,1], s=70, label="ORCA", alpha=0.85)
    ax.scatter(coords_ref[:,0], coords_ref[:,1], s=70, marker="x", label="NH₃ (ref)")
    ax.set_title("Proyección XY")
    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)")
    ax.grid(alpha=0.3); ax.legend()
    st.subheader("Gráfica XY")
    st.pyplot(fig)

    # Tabla de métricas
    fila_orca = _fila_metricas("ORCA", elems_o, coords_o)
    fila_nh3  = _fila_metricas("NH₃", elems_ref, coords_ref)
    df_m = pd.DataFrame([fila_orca, fila_nh3], columns=[
        "Molécula","Átomos","Masa (u)","Elementos","Centro (Å)","Caja Δx,Δy,Δz (Å)",
        "R_g (Å)","Distancia media (Å)","Distancia máx (Å)","⟨N–H⟩ (Å)","⟨C–H⟩ (Å)"
    ])
    st.subheader("Tabla comparativa")
    st.dataframe(df_m, use_container_width=True)
    return df_m

#entrada para la optimización de las geometrías moleculares
def comparar_moleculas_orca_vs_agua(ruta_paso_3):
    simbolos, coords = [], []
    leer = False
    with open(ruta_paso_3, "r") as f:
        for linea in f:
            linea = linea.strip()
            if linea.lower().startswith("* xyz"):
                leer = True
                continue
            if leer:
                if linea.startswith("*"):
                    break
                partes = linea.split()
                if len(partes) == 4:
                    simbolos.append(partes[0])
                    coords.append([float(x) for x in partes[1:]])

    # Agua de referencia
    simbolos_h2o = ["O", "H", "H"]
    coords_h2o = [
        [0.000, 0.000, 0.000],   # O
        [0.958, 0.000, 0.000],   # H
        [-0.239, 0.927, 0.000]   # H
    ]

    # ---- Gráfico 3D ----
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Molécula Orca
    xs1, ys1, zs1 = zip(*coords)
    ax.scatter(xs1, ys1, zs1, c="blue", label="Molécula Orca", s=60)
    for sym, (x, y, z) in zip(simbolos, coords):
        ax.text(x, y, z, sym, fontsize=8, color="blue")

    # Agua
    xs2, ys2, zs2 = zip(*coords_h2o)
    ax.scatter(xs2, ys2, zs2, c="red", label="H2O", s=60)
    for sym, (x, y, z) in zip(simbolos_h2o, coords_h2o):
        ax.text(x, y, z, sym, fontsize=8, color="red")

    ax.set_title("Comparación: Molécula Orca vs Agua")
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.legend()

    st.pyplot(fig)

    # ---- Tabla comparativa ----
    # Ajustamos el tamaño para que ambas listas coincidan
    max_len = max(len(simbolos), len(simbolos_h2o))
    simbolos += ["-"] * (max_len - len(simbolos))
    coords += [[None, None, None]] * (max_len - len(coords))
    simbolos_h2o += ["-"] * (max_len - len(simbolos_h2o))
    coords_h2o += [[None, None, None]] * (max_len - len(coords_h2o))

    data = []
    for s_orca, c_orca, s_h2o, c_h2o in zip(simbolos, coords, simbolos_h2o, coords_h2o):
        data.append({
            "Átomo Orca": s_orca,
            "X_Orca": c_orca[0], "Y_Orca": c_orca[1], "Z_Orca": c_orca[2],
            "Átomo H2O": s_h2o,
            "X_H2O": c_h2o[0], "Y_H2O": c_h2o[1], "Z_H2O": c_h2o[2],
        })

    df = pd.DataFrame(data)
    st.subheader("Tabla comparativa de coordenadas")
    st.dataframe(df, use_container_width=True)

# comparar IR y Raman (por modos) by Francis  
def comparar_ir_raman_vs_nh3(ruta_paso_3):
    # -----------------------
    # Lectura de datos ORCA
    # -----------------------
    if not os.path.exists(ruta_paso_3):
        st.error(f"No se encontró el archivo:\n{ruta_paso_3}")
        return

    ir_data, raman_data = [], []
    seccion = None
    with open(ruta_paso_3, "r") as f:
        for linea in f:
            linea = linea.strip()
            if not linea or linea.startswith("#"):
                continue
            if linea.startswith("IR SPECTRUM"):
                seccion = "IR";  continue
            elif linea.startswith("RAMAN SPECTRUM"):
                seccion = "RAMAN";  continue
            elif linea[0].isdigit():  # datos numéricos
                partes = (linea.replace(":", "")
                                .replace("(", "")
                                .replace(")", "")
                                .split())
                if seccion == "IR" and len(partes) >= 8:
                    modo, freq, eps, intensidad, t2, tx, ty, tz = partes[:8]
                    ir_data.append({
                        "Mode": int(modo),
                        "Freq (cm^-1)": float(freq),
                        "Int (km/mol)": float(intensidad)
                    })
                elif seccion == "RAMAN" and len(partes) >= 4:
                    modo, freq, actividad, depol = partes[:4]
                    raman_data.append({
                        "Mode": int(modo),
                        "Freq (cm^-1)": float(freq),
                        "Actividad": float(actividad)
                    })

    df_ir = pd.DataFrame(ir_data)
    df_raman = pd.DataFrame(raman_data)

    # -----------------------
    # Datos de referencia NH3
    # (aprox. fundam.: ν1~3336, ν2~950, ν3~3444, ν4~1627 cm-1)
    # Intensidades/actividades son ilustrativas (ajústalas si tienes valores de tu referencia)
    # -----------------------
    df_ir_nh3 = pd.DataFrame({
        "Mode": [1, 2, 3, 4],
        "Freq (cm^-1)": [950, 1627, 3336, 3444],
        "Int (km/mol)": [0.4, 0.8, 1.2, 1.0]
    })

    df_raman_nh3 = pd.DataFrame({
        "Mode": [1, 2, 3, 4],
        "Freq (cm^-1)": [950, 1627, 3336, 3444],
        "Actividad": [0.9, 0.7, 0.5, 0.8]
    })

    # -----------------------
    # Graficar IR por modos
    # -----------------------
    fig_ir, ax_ir = plt.subplots(figsize=(7, 4))
    if not df_ir.empty and "Mode" in df_ir and "Int (km/mol)" in df_ir:
        ax_ir.plot(df_ir["Mode"], df_ir["Int (km/mol)"], "-o", label="Molécula ORCA")
    else:
        st.warning("IR: No se encontraron columnas esperadas en ORCA. Se requieren 'Mode' e 'Int (km/mol)'.")
    ax_ir.plot(df_ir_nh3["Mode"], df_ir_nh3["Int (km/mol)"], "-s", label="NH₃ (ref)")
    ax_ir.set_title("Comparación Espectro IR (por modo)")
    ax_ir.set_xlabel("Modo vibracional")
    ax_ir.set_ylabel("Intensidad (km/mol)")
    ax_ir.legend()

    # -----------------------
    # Graficar Raman por modos
    # -----------------------
    fig_raman, ax_raman = plt.subplots(figsize=(7, 4))
    if not df_raman.empty and "Mode" in df_raman and "Actividad" in df_raman:
        ax_raman.plot(df_raman["Mode"], df_raman["Actividad"], "-o", label="Molécula ORCA")
    else:
        st.warning("Raman: No se encontraron columnas esperadas en ORCA. Se requieren 'Mode' y 'Actividad'.")
    ax_raman.plot(df_raman_nh3["Mode"], df_raman_nh3["Actividad"], "-s", label="NH₃ (ref)")
    ax_raman.set_title("Comparación Espectro Raman (por modo)")
    ax_raman.set_xlabel("Modo vibracional")
    ax_raman.set_ylabel("Actividad")
    ax_raman.legend()

    # -----------------------
    # Mostrar en Streamlit: Tablas y Gráficas
    # -----------------------
    # 1) Comparación por Mode (igual que tu flujo original)
    st.subheader("Tabla IR (Molécula vs NH₃) – por Mode")
    if not df_ir.empty and "Mode" in df_ir and "Int (km/mol)" in df_ir:
        st.dataframe(
            df_ir.merge(df_ir_nh3[["Mode", "Int (km/mol)"]],
                        on="Mode", how="outer",
                        suffixes=("_Orca", "_NH3"))
        )
    else:
        st.info("No se puede comparar IR por 'Mode' (faltan columnas en ORCA).")

    st.subheader("Comparación Espectro IR")
    st.pyplot(fig_ir)

    st.subheader("Tabla Raman (Molécula vs NH₃) – por Mode")
    if not df_raman.empty and "Mode" in df_raman and "Actividad" in df_raman:
        st.dataframe(
            df_raman.merge(df_raman_nh3[["Mode", "Actividad"]],
                           on="Mode", how="outer",
                           suffixes=("_Orca", "_NH3"))
        )
    else:
        st.info("No se puede comparar Raman por 'Mode' (faltan columnas en ORCA).")

    st.subheader("Comparación Espectro Raman")
    st.pyplot(fig_raman)

    # -----------------------
    # PLUS: Comparación por frecuencia más cercana (útil si los Mode no coinciden)
    # -----------------------
    def join_por_frecuencia(df_a, df_b, col_y_a, col_y_b, tol=50.0):
        """
        Para cada frecuencia de df_a, busca la más cercana en df_b.
        Si la diferencia <= tol (cm^-1), la empareja; si no, queda NaN.
        """
        if df_a.empty or df_b.empty:
            return pd.DataFrame()
        out = []
        for _, r in df_a.iterrows():
            fa = r["Freq (cm^-1)"]
            idx_min = (df_b["Freq (cm^-1)"] - fa).abs().idxmin()
            fb = df_b.loc[idx_min, "Freq (cm^-1)"]
            diff = abs(fb - fa)
            ya = r.get(col_y_a, None)
            yb = df_b.loc[idx_min, col_y_b]
            out.append({"Freq_ORCA": fa, col_y_a: ya,
                        "Freq_NH3": fb, col_y_b: yb, "Δ|cm^-1|": diff})
        df_out = pd.DataFrame(out)
        return df_out[df_out["Δ|cm^-1|"] <= tol].reset_index(drop=True)

    st.subheader("IR – Emparejamiento por frecuencia cercana (±50 cm⁻¹)")
    if not df_ir.empty:
        df_ir_freq = join_por_frecuencia(df_ir, df_ir_nh3, "Int (km/mol)", "Int (km/mol)", tol=50.0)
        if df_ir_freq.empty:
            st.info("No hubo coincidencias IR dentro de ±50 cm⁻¹.")
        else:
            st.dataframe(df_ir_freq, use_container_width=True)
    else:
        st.info("No hay datos IR de ORCA para emparejar por frecuencia.")

    st.subheader("Raman – Emparejamiento por frecuencia cercana (±50 cm⁻¹)")
    if not df_raman.empty:
        df_ra_freq = join_por_frecuencia(df_raman, df_raman_nh3, "Actividad", "Actividad", tol=50.0)
        if df_ra_freq.empty:
            st.info("No hubo coincidencias Raman dentro de ±50 cm⁻¹.")
        else:
            st.dataframe(df_ra_freq, use_container_width=True)
    else:
        st.info("No hay datos Raman de ORCA para emparejar por frecuencia.")

# ----- Comparación S4 (RMN) vs NH3 -----
def comparar_rmn_s4_vs_nh3(ruta: str) -> pd.DataFrame:
    """
    Mismo estilo: lee ORCA_input_RMN_S4 (bloque * xyz ... *),
    muestra tablas/gráfica y retorna tabla de métricas.
    """
    st.header("RMN S4 · Comparación geometría ORCA vs NH₃")

    if not os.path.exists(ruta):
        st.error(f"No se encontró el archivo:\n{ruta}")
        return pd.DataFrame()

    # Parsear
    elems_o, coords_o = _parse_orca_xyz(ruta)
    elems_ref, coords_ref = _nh3_ref()

    # Tablas de coordenadas
    st.subheader("Coordenadas ORCA (S4)")
    st.dataframe(_df_coordenadas("ORCA (S4)", elems_o, coords_o), use_container_width=True)

    st.subheader("Coordenadas NH₃ (referencia)")
    st.dataframe(_df_coordenadas("NH₃", elems_ref, coords_ref), use_container_width=True)

    # Gráfica
    fig, ax = plt.subplots(figsize=(7,4))
    if coords_o.size>0:
        ax.scatter(coords_o[:,0], coords_o[:,1], s=70, c="tab:green", label="ORCA (S4)", alpha=0.85)
    ax.scatter(coords_ref[:,0], coords_ref[:,1], s=70, c="tab:purple", marker="x", label="NH₃ (ref)")
    ax.set_title("Proyección XY")
    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)")
    ax.grid(alpha=0.3); ax.legend()
    st.subheader("Gráfica XY")
    st.pyplot(fig)

    # Tabla de métricas
    fila_orca = _fila_metricas("ORCA (S4)", elems_o, coords_o)
    fila_nh3  = _fila_metricas("NH₃", elems_ref, coords_ref)
    df_m = pd.DataFrame([fila_orca, fila_nh3], columns=[
        "Molécula","Átomos","Masa (u)","Elementos","Centro (Å)","Caja Δx,Δy,Δz (Å)",
        "R_g (Å)","Distancia media (Å)","Distancia máx (Å)","⟨N–H⟩ (Å)","⟨C–H⟩ (Å)"
    ])
    st.subheader("Tabla comparativa")
    st.dataframe(df_m, use_container_width=True)
    return df_m

# ---------- MAIN ----------
def main():
    # Menú de navegación
    option = st.sidebar.radio("Selecciona una sección", [
    "Molecula 3D", "Molécula 2D", "Conjunto de moléculas", "Contenedor de moleculas", 
    "Espectro IR", "Trabajo de adhesión porcentual", "Molécula teórica","Raman", "Comparacion Moleculas",
    "Raman vs molecula", "Estimacion de desplazamiento"
    ])

    if option == "Molecula 3D":
        fig = dibujar_nh3()
        st.title("Molécula de NH₃ en 3D")
        st.write("Visualización interactiva usando Plotly y Streamlit")
        st.plotly_chart(fig, use_container_width=True)
        
    elif option == "Molécula 2D":
        st.title("Molécula")
        fig = dibujar_nh3_2d()
        st.pyplot(fig)
        
    elif option == "Conjunto de moléculas":
        st.title("Conjunto de moléculas")
        fig = dibujar_conjunto_nh3()
        st.plotly_chart(fig, use_container_width=True)
        
    elif option == "Contenedor de moleculas":
        st.title("Contenedor de moléculas")
        fig = contenedor_molecular()
        st.plotly_chart(fig, use_container_width=True)

    elif option == "Espectro IR":
        st.title("Espectro IR")
        fig = dibujar_ir()
        st.pyplot(fig)
        
    elif option == "Trabajo de adhesión porcentual":
        st.title("Trabajo de adhesión porcentual")
        fig = graficar_trabajo_adhesion()
        st.pyplot(fig)

    elif option == "Molécula teórica":
        st.title("Molécula teórica")
        mostrar_rdf()

    elif option == "Comparacion Moleculas":
        st.title("Comparar Moleculas ORCA vs NH3")
        ruta = "modelos/paso_2.txt"
        comparar_moleculas_orca_vs_nh3(ruta)

    elif option == "Raman":
        st.title("Graficos IR y Raman")
        ruta_paso_3 = "modelos/paso_3.txt"
        mostrar_ir_raman(ruta_paso_3)

    elif option == "Raman vs molecula":
        st.title("Graficos IR y Raman comparado con molecula NH3")
        ruta = "modelos/paso_3.txt" 
        comparar_ir_raman_vs_nh3(ruta)
    
    elif option == "Estimacion de desplazamiento":
        st.title("Estimacion De Los Desplazamientos Químicos de H1 y C13")
        ruta = "modelos/paso_4.txt"
        comparar_rmn_s4_vs_nh3(ruta)


if __name__ == "__main__":
    main()
