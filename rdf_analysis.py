# rdf_analysis.py
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf
import matplotlib.pyplot as plt
from pathlib import Path
import streamlit as st
import pandas as pd

def generate_trajectory(molecule_name, n_frames=100, n_molecules=50, box_size=20.0):
    """
    Genera una trayectoria de prueba con varias moléculas.
    """
    try:
        from molecular_config_manager import MolecularConfigManager
        config_manager = MolecularConfigManager()
        elements, coordinates, _ = config_manager.read_xyz_file(molecule_name)

        if len(elements) == 0:
            st.error("No se pudieron obtener las coordenadas de la molécula")
            return None

        temp_dir = Path("temp_trajectories")
        temp_dir.mkdir(exist_ok=True)
        traj_file = temp_dir / f"{molecule_name}_traj.xyz"

        # Construir varias copias aleatorias en la caja
        with open(traj_file, 'w') as f:
            for frame in range(n_frames):
                n_atoms_total = len(elements) * n_molecules
                f.write(f'{n_atoms_total}\n')
                f.write(f'Frame: {frame} - {molecule_name}\n')

                for m in range(n_molecules):
                    # offset aleatorio dentro de la caja
                    offset = np.random.uniform(0, box_size, 3)
                    for element, coord in zip(elements, coordinates):
                        noise = np.random.normal(0, 0.05, 3)  # vibración pequeña
                        new_coord = coord + offset + noise
                        f.write(f'{element} {new_coord[0]:.6f} {new_coord[1]:.6f} {new_coord[2]:.6f}\n')

        st.success(f"Trayectoria generada con {n_molecules} moléculas: {traj_file}")
        return str(traj_file)

    except Exception as e:
        st.error(f"Error generando trayectoria: {e}")
        return None

def calculate_rdf(traj_file, element_pairs, range=(0.5, 5.0), nbins=100, box_size=20.0):
    """
    Calcula la función de distribución radial para pares de elementos específicos.
    Excluye distancias intramoleculares.
    """
    try:
        u = mda.Universe(traj_file, format='XYZ')
        u.dimensions = np.array([box_size, box_size, box_size, 90, 90, 90])

        results = {}

        for element1, element2 in element_pairs:
            atoms1 = u.select_atoms(f'element {element1}')
            atoms2 = u.select_atoms(f'element {element2}')

            if len(atoms1) == 0 or len(atoms2) == 0:
                st.warning(f"No se encontraron átomos para el par {element1}-{element2}")
                continue

            # exclusion_block=(1,1) evita contar átomos del mismo residuo/molécula
            rdf_calc = rdf.InterRDF(atoms1, atoms2, range=range, nbins=nbins, exclusion_block=(1, 1))
            rdf_calc.run()

            results[f"{element1}-{element2}"] = {
                'bins': rdf_calc.results.bins,
                'rdf': rdf_calc.results.rdf,
                'counts': rdf_calc.results.count
            }

        return results

    except Exception as e:
        st.error(f"Error calculando RDF: {e}")
        return None

def plot_rdf(results, molecule_name):
    """
    Genera gráficos de las funciones de distribución radial.
    """
    if not results:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, len(results)))

    for i, (pair, data) in enumerate(results.items()):
        ax.plot(data['bins'], data['rdf'],
                color=colors[i], linewidth=2, label=pair)

    ax.set_xlabel('Distancia (Å)')
    ax.set_ylabel('g(r)')
    ax.set_title(f'Función de Distribución Radial - {molecule_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    return fig

def get_suggested_params(molecule_name, n_molecules):
    """
    Devuelve parámetros sugeridos (frames, box_size, r_min, r_max, n_bins)
    según la molécula seleccionada y número de moléculas.
    """
    # densidades aproximadas en moléculas/Å³
    densities = {
        "water": 0.033456,  # 1 g/cm3
        "methanol": 0.026,  # 0.79 g/cm3
        "ethanol": 0.021,   # 0.79 g/cm3 aprox.
    }

    density = densities.get(molecule_name.lower(), 0.033456)
    box_size = (n_molecules / density) ** (1/3)

    return {
        "n_frames": 100,
        "box_size": box_size,
        "r_min": 0.5,
        "r_max": 6.0,
        "n_bins": 150
    }


def analyze_molecule_rdf(molecule_name):
    
    """
    Función principal para analizar la RDF de una molécula.
    """
    st.header(f"📊 Análisis de Función de Distribución Radial - {molecule_name}")

    n_molecules = st.slider("Número de moléculas", 1, 200, 50)

    # Obtener sugerencias
    suggested = get_suggested_params(molecule_name, n_molecules)

    st.info(
        f"💡 Sugerencia para {molecule_name}: "
        f"{suggested['n_frames']} frames, caja {suggested['box_size']:.2f} Å, "
        f"distancia min {suggested['r_min']} Å, max {suggested['r_max']} Å, "
        f"{suggested['n_bins']} bins."
    )

    col1, col2 = st.columns(2)

    with col1:
        n_frames = st.slider("Número de frames", 10, 500, suggested["n_frames"])
        box_size = st.slider("Tamaño de caja (Å)", 10.0, 100.0, float(suggested["box_size"]))

    with col2:
        r_min = st.slider("Distancia mínima (Å)", 0.1, 2.0, suggested["r_min"])
        r_max = st.slider("Distancia máxima (Å)", 3.0, 10.0, suggested["r_max"])
        n_bins = st.slider("Número de bins", 50, 300, suggested["n_bins"])


    # Obtener elementos de la molécula
    from molecular_config_manager import MolecularConfigManager
    config_manager = MolecularConfigManager()
    elements, _, _ = config_manager.read_xyz_file(molecule_name)

    if not elements:
        st.error("No se pudieron obtener los elementos de la molécula")
        return

    unique_elements = list(set(elements))

    # Seleccionar pares de elementos para analizar
    st.subheader("Seleccionar pares de elementos")
    selected_pairs = []

    for i, elem1 in enumerate(unique_elements):
        for j, elem2 in enumerate(unique_elements[i:], i):
            if st.checkbox(f"{elem1}-{elem2}", value=(i == j), key=f"pair_{elem1}_{elem2}"):
                selected_pairs.append((elem1, elem2))

    if not selected_pairs:
        st.warning("Selecciona al menos un par de elementos")
        return

    # Ejecutar análisis
    if st.button("Calcular RDF", type="primary"):
        with st.spinner("Generando trayectoria y calculando RDF..."):
            # Generar trayectoria con varias moléculas
            traj_file = generate_trajectory(
                molecule_name,
                n_frames=n_frames,
                n_molecules=n_molecules,
                box_size=box_size
            )

            if traj_file:
                # Calcular RDF
                results = calculate_rdf(
                    traj_file,
                    selected_pairs,
                    range=(r_min, r_max),
                    nbins=n_bins,
                    box_size=box_size
                )

                if results:
                    # Mostrar resultados
                    st.subheader("Resultados")

                    # Gráfico
                    fig = plot_rdf(results, molecule_name)
                    if fig:
                        st.pyplot(fig)

                    # Tabla de valores
                    st.subheader("Valores de RDF")

                    for pair, data in results.items():
                        with st.expander(f"Valores para {pair}"):
                            df = pd.DataFrame({
                                'Distancia (Å)': data['bins'],
                                'g(r)': data['rdf'],
                                'Conteos': data['counts']
                            })
                            st.dataframe(df, use_container_width=True)

                    # Información estadística
                    st.subheader("Estadísticas")

                    cols = st.columns(len(results))

                    for i, (pair, data) in enumerate(results.items()):
                        with cols[i]:
                            st.metric(f"{pair} - Máximo", f"{np.max(data['rdf']):.3f}")
                            st.metric(f"{pair} - Distancia pico",
                                    f"{data['bins'][np.argmax(data['rdf'])]:.2f} Å")

                    # Botón para descargar datos
                    csv_data = "Pair,Distance (Å),g(r),Counts\n"

                    for pair, data in results.items():
                        for dist, g_val, count in zip(data['bins'], data['rdf'], data['counts']):
                            csv_data += f"{pair},{dist:.4f},{g_val:.4f},{count}\n"

                    st.download_button(
                        label="Descargar datos CSV",
                        data=csv_data,
                        file_name=f"rdf_{molecule_name}.csv",
                        mime="text/csv"
                    )
                else:
                    st.error("Error calculando RDF")
