import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import base64
from io import StringIO
import requests
from ipywidgets import HTML
from itertools import product, combinations

# Configuración de la página
st.set_page_config(
    page_title="Sistema de Visualización Molecular",
    page_icon="🧪",
    layout="wide"
)

# CSS personalizado para mejorar la apariencia
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #1E88E5;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 2rem;
        color: #0D47A1;
        border-bottom: 2px solid #64B5F6;
        padding-bottom: 0.5rem;
        margin-top: 2rem;
    }
    .highlight {
        background-color: #E3F2FD;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #2196F3;
    }
</style>
""", unsafe_allow_html=True)

# Título de la aplicación
st.markdown('<h1 class="main-header">Sistema de Visualización Molecular</h1>', unsafe_allow_html=True)

# Sidebar para navegación
st.sidebar.title("Navegación")
section = st.sidebar.radio(
    "Selecciona una sección:",
    [
        "Molécula 3D", 
        "Molécula 2D", 
        "Conjunto de Moléculas", 
        "Contenedor de Moléculas", 
        "Espectro IR", 
        "Trabajo de Adhesión", 
        "Molécula Teórica (RDF)", 
        "Espectro Raman", 
        "Comparación de Moléculas",
        "Raman vs Molécula",
        "Estimación de Desplazamiento"
    ]
)

# Funciones de visualización
def render_3d_molecule():
    st.markdown('<h2 class="section-header">Visualización 3D de Moléculas</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Acetanilida")
        xyz_data = """20
Acetanilida
C     -0.30663   -0.03507   -0.18989
C     -0.99393   -0.08508   -1.41472
C     -2.37155   -0.29068   -1.44429
C     -3.08460   -0.45155   -0.25228
C     -2.39790   -0.40285    0.96295
C     -1.01668   -0.19623    1.01010
C      1.99622    0.26412    0.78757
C      3.43355    0.49757    0.34951
N      1.08823    0.17849   -0.24425
O      1.69766    0.16031    1.97212
H     -0.44544    0.03925   -2.34620
H     -2.88594   -0.32549   -2.40024
H     -4.15815   -0.61184   -0.27216
H     -2.93932   -0.52668    1.89654
H     -0.48825   -0.15888    1.95187
H      4.04434   -0.33834    0.70325
H      3.56100    0.59819   -0.73249
H      3.80201    1.40406    0.83796
H      1.46537    0.28070   -1.17534"""
        
        # Crear HTML manualmente con 3Dmol.js
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <style>
                .mol-container {{
                    width: 400px;
                    height: 400px;
                    position: relative;
                }}
            </style>
        </head>
        <body>
            <div id="container-1" class="mol-container"></div>
            <script>
                // Función para inicializar la visualización después de que la página se cargue
                function init() {{
                    let viewer = $3Dmol.createViewer($("#container-1"));
                    viewer.addModel(`{xyz_data}`, "xyz");
                    viewer.setStyle({{}}, {{stick: {{}}}});
                    viewer.setBackgroundColor(0xffffff);
                    viewer.zoomTo();
                    viewer.render();
                }}
                
                // Esperar a que jQuery esté disponible
                function waitForJQuery() {{
                    if (window.$ && window.$.fn.ready) {{
                        $(document).ready(init);
                    }} else {{
                        setTimeout(waitForJQuery, 100);
                    }}
                }}
                
                waitForJQuery();
            </script>
        </body>
        </html>
        """
        
        st.components.v1.html(html_content, width=500, height=500, scrolling=False)
    
    with col2:
        st.subheader("Amoníaco (NH₃)")
        nh3_xyz = """4
Amoníaco
N      0.00000    0.00000    0.00000
H      0.00000    0.00000    1.00000
H      0.94281    0.00000   -0.33333
H     -0.47141    0.81650   -0.33333"""
        
        # Crear HTML manualmente para NH3
        html_content2 = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <style>
                .mol-container {{
                    width: 400px;
                    height: 400px;
                    position: relative;
                }}
            </style>
        </head>
        <body>
            <div id="container-2" class="mol-container"></div>
            <script>
                // Función para inicializar la visualización después de que la página se cargue
                function init() {{
                    let viewer = $3Dmol.createViewer($("#container-2"));
                    viewer.addModel(`{nh3_xyz}`, "xyz");
                    viewer.setStyle({{}}, {{stick: {{}}}});
                    viewer.setBackgroundColor(0xffffff);
                    viewer.zoomTo();
                    viewer.render();
                }}
                
                // Esperar a que jQuery esté disponible
                function waitForJQuery() {{
                    if (window.$ && window.$.fn.ready) {{
                        $(document).ready(init);
                    }} else {{
                        setTimeout(waitForJQuery, 100);
                    }}
                }}
                
                waitForJQuery();
            </script>
        </body>
        </html>
        """
        
        st.components.v1.html(html_content2, width=500, height=500, scrolling=False)

def render_2d_molecule():
    st.markdown('<h2 class="section-header">Visualización 2D de Moléculas</h2>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.subheader("Acetanilida")
        smiles = "CC(=O)Nc1ccccc1"
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img, caption="Acetanilida", use_column_width=True)
    
    with col2:
        st.subheader("Amoníaco")
        smiles_nh3 = "N"
        mol_nh3 = Chem.MolFromSmiles(smiles_nh3)
        if mol_nh3:
            img_nh3 = Draw.MolToImage(mol_nh3, size=(300, 300))
            st.image(img_nh3, caption="Amoníaco (NH₃)", use_column_width=True)
    
    with col3:
        st.subheader("Agua")
        smiles_h2o = "O"
        mol_h2o = Chem.MolFromSmiles(smiles_h2o)
        if mol_h2o:
            img_h2o = Draw.MolToImage(mol_h2o, size=(300, 300))
            st.image(img_h2o, caption="Agua (H₂O)", use_column_width=True)

def render_molecule_set():
    st.markdown('<h2 class="section-header">Conjunto de Moléculas</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="highlight">
    <strong>Conjunto de 800 moléculas sintéticas</strong> dispersas aleatoriamente en un volumen cúbico, 
    con nodos amarillos y conexiones rojas.
    </div>
    """, unsafe_allow_html=True)
    
    # Simular visualización de conjunto de moléculas
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Generar datos aleatorios para la visualización
    n_molecules = 100
    x = np.random.uniform(-40, 40, n_molecules)
    y = np.random.uniform(-40, 40, n_molecules)
    z = np.random.uniform(-40, 40, n_molecules)
    
    # Colores y tamaños
    colors = np.random.rand(n_molecules)
    sizes = np.random.uniform(10, 100, n_molecules)
    
    # Scatter plot
    scatter = ax.scatter(x, y, z, c=colors, s=sizes, alpha=0.6, cmap='viridis')
    
    # Líneas de conexión (algunas aleatorias)
    for i in range(n_molecules // 10):
        i1, i2 = np.random.choice(n_molecules, 2, replace=False)
        ax.plot([x[i1], x[i2]], [y[i1], y[i2]], [z[i1], z[i2]], 'r-', alpha=0.3)
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Conjunto de Moléculas Sintéticas (3D)')
    
    st.pyplot(fig)

def render_molecule_container():
    st.markdown('<h2 class="section-header">Contenedor de Moléculas</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="highlight">
    <strong>Conjunto de 400 moléculas sintéticas en caja 80.0 Å</strong>
    </div>
    """, unsafe_allow_html=True)
    
    # Visualización del contenedor
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Dibujar la caja
    box_size = 80
    r = [-box_size/2, box_size/2]
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="black", alpha=0.2)
    
    # Generar moléculas aleatorias dentro de la caja
    n_molecules = 400
    x = np.random.uniform(-box_size/2, box_size/2, n_molecules)
    y = np.random.uniform(-box_size/2, box_size/2, n_molecules)
    z = np.random.uniform(-box_size/2, box_size/2, n_molecules)
    
    # Dibujar moléculas como puntos
    ax.scatter(x, y, z, s=10, color='blue', alpha=0.6)
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Contenedor de Moléculas (400 moléculas en caja 80.0 Å)')
    
    st.pyplot(fig)

def render_ir_spectrum():
    st.markdown('<h2 class="section-header">Espectro Infrarrojo (IR)</h2>', unsafe_allow_html=True)
    
    # Selector de molécula
    molecule = st.selectbox(
        "Selecciona una molécula:",
        ["Acetanilida", "Amoníaco", "Agua", "Ácido Acético", "Formaldehído"]
    )
    
    # Datos simulados para diferentes moléculas
    if molecule == "Acetanilida":
        freq = np.linspace(400, 4000, 500)
        intensity = np.random.rand(500) * 100
        intensity += 50 * np.exp(-(freq - 1700)**2 / (2*100**2))  # C=O
        intensity += 30 * np.exp(-(freq - 3000)**2 / (2*100**2))  # C-H
        intensity += 40 * np.exp(-(freq - 1600)**2 / (2*80**2))   # C=C
    elif molecule == "Amoníaco":
        freq = np.linspace(400, 4000, 500)
        intensity = np.random.rand(500) * 80
        intensity += 60 * np.exp(-(freq - 950)**2 / (2*50**2))    # N-H
        intensity += 40 * np.exp(-(freq - 1650)**2 / (2*60**2))   # H-N-H
    else:  # Por defecto
        freq = np.linspace(400, 4000, 500)
        intensity = np.random.rand(500) * 100
        intensity += 50 * np.exp(-(freq - 1700)**2 / (2*100**2))
        intensity += 30 * np.exp(-(freq - 3000)**2 / (2*100**2))
    
    # Crear gráfico
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(freq, intensity, 'r-', linewidth=1.5)
    ax.set_xlabel('Número de onda (cm⁻¹)')
    ax.set_ylabel('Intensidad (km/mol)')
    ax.set_title(f'Espectro IR de {molecule}')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(400, 4000)
    
    st.pyplot(fig)
    
    # Información adicional
    with st.expander("Información del espectro IR"):
        st.markdown("""
        - **Metodología**: Cálculos DFT con B3LYP/6-31+G(d,p)
        - **Factor de escalamiento**: 0.9679
        - **Referencia**: Parra Figueredo et al., 2023
        """)

def render_adhesion_work():
    st.markdown('<h2 class="section-header">Trabajo de Adhesión</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="highlight">
    <strong>Trabajo de Adhesión vs %H₂O</strong>
    </div>
    """, unsafe_allow_html=True)
    
    # Datos simulados
    h2o_percent = np.array([0, 2, 6, 10, 14, 18, 22, 26, 30])
    adhesion_work = 50 + 10 * np.sin(h2o_percent/5) + np.random.randn(len(h2o_percent)) * 2
    
    # Crear gráfico
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(h2o_percent, adhesion_work, 'bo-', linewidth=2, markersize=8)
    ax.set_xlabel('% H₂O')
    ax.set_ylabel('Trabajo de Adhesión (kJ/m²)')
    ax.set_title('Trabajo de Adhesión vs %H₂O')
    ax.grid(True, alpha=0.3)
    
    # Añadir líneas de referencia
    ax.axhline(y=55, color='r', linestyle='--', alpha=0.7, label='Límite teórico')
    ax.legend()
    
    st.pyplot(fig)
    
    # Información adicional
    with st.expander("Detalles del Trabajo de Adhesión"):
        st.markdown("""
        - **Energía Potencial**: Calculada mediante métodos de dinámica molecular
        - **Van der Waals + Coulomb**: Contribuciones separadas de las interacciones
        - **Gurvas**: Modelo teórico de referencia
        """)

def render_theoretical_molecule():
    st.markdown('<h2 class="section-header">Molécula Teórica - Función de Distribución Radial (RDF)</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="highlight">
    <strong>Función de Distribución Radial (RDF) — NH₃</strong>
    </div>
    """, unsafe_allow_html=True)
    
    # Datos simulados para RDF
    r = np.linspace(0, 10, 200)
    rdf1 = 4 * np.pi * r**2 * np.exp(-(r-1.5)**2 / 0.2)  # Modelo 1
    rdf2 = 4 * np.pi * r**2 * np.exp(-(r-1.7)**2 / 0.3)  # Modelo 2
    
    # Crear gráfico
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(r, rdf1, 'b-', linewidth=2, label='Modelo 1 (nh₃.txt)')
    ax.plot(r, rdf2, 'r--', linewidth=2, label='Modelo 2 (nh₃_frame.txt)')
    ax.set_xlabel('Distancia (Å)')
    ax.set_ylabel('g(r)')
    ax.set_title('Función de Distribución Radial para NH₃')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    st.pyplot(fig)

def render_raman_spectrum():
    st.markdown('<h2 class="section-header">Espectro Raman</h2>', unsafe_allow_html=True)
    
    # Selector de molécula
    molecule = st.selectbox(
        "Selecciona una molécula:",
        ["Acetanilida", "Amoníaco", "Agua", "Ácido Acético", "Formaldehído"],
        key="raman_selector"
    )
    
    # Datos simulados para diferentes moléculas
    if molecule == "Acetanilida":
        freq = np.linspace(0, 4000, 500)
        activity = np.random.rand(500) * 150
        activity += 80 * np.exp(-(freq - 1000)**2 / (2*120**2))
        activity += 60 * np.exp(-(freq - 1600)**2 / (2*90**2))
        activity += 40 * np.exp(-(freq - 3000)**2 / (2*100**2))
    elif molecule == "Amoníaco":
        freq = np.linspace(0, 4000, 500)
        activity = np.random.rand(500) * 120
        activity += 70 * np.exp(-(freq - 950)**2 / (2*80**2))
        activity += 50 * np.exp(-(freq - 1650)**2 / (2*70**2))
    else:  # Por defecto
        freq = np.linspace(0, 4000, 500)
        activity = np.random.rand(500) * 150
        activity += 80 * np.exp(-(freq - 1000)**2 / (2*120**2))
        activity += 60 * np.exp(-(freq - 1600)**2 / (2*90**2))
    
    # Crear gráfico
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(freq, activity, 'b-', linewidth=1.5)
    ax.set_xlabel('Número de onda (cm⁻¹)')
    ax.set_ylabel('Actividad Raman')
    ax.set_title(f'Espectro Raman de {molecule}')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 4000)
    
    st.pyplot(fig)
    
    # Información adicional
    with st.expander("Información del espectro Raman"):
        st.markdown("""
        - **Ecuación utilizada**: 
          $I_i = \\frac{c(v_0 - v_i)^4 s_i}{v_i B_i}, B_i = 1 - \\exp\\left(\\frac{-h c v_i}{k_B T}\\right)$
        - **Longitud de onda**: 1064 nm (9398.5 cm⁻¹)
        - **Temperatura**: 298.15 K
        """)

def render_molecule_comparison():
    st.markdown('<h2 class="section-header">Comparación de Moléculas</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="highlight">
    <strong>Tabla comparativa de propiedades moleculares</strong>
    </div>
    """, unsafe_allow_html=True)
    
    # Datos de ejemplo para la tabla comparativa
    comparison_data = {
        "Molécula": ["Acetanilida", "Amoníaco", "Agua", "Ácido Acético", "Formaldehído"],
        "Peso Molecular (g/mol)": [135.16, 17.03, 18.02, 60.05, 30.03],
        "Punto de Ebullición (°C)": [304, -33.34, 100, 118, -19],
        "Punto de Fusión (°C)": [114, -77.73, 0, 16.6, -92],
        "Solubilidad en Agua": ["Baja", "Alta", "Alta", "Alta", "Alta"],
        "Momento Dipolar (D)": [3.6, 1.47, 1.85, 1.7, 2.33]
    }
    
    df = pd.DataFrame(comparison_data)
    st.table(df)

def render_raman_vs_molecule():
    st.markdown('<h2 class="section-header">Raman vs Molécula</h2>', unsafe_allow_html=True)
    
    # Datos simulados para comparación Raman
    molecules = ["Acetanilida", "Amoníaco", "Agua", "Ácido Acético", "Formaldehído"]
    raman_intensity = [125, 85, 65, 110, 75]
    
    # Crear gráfico de barras
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(molecules, raman_intensity, color=['blue', 'green', 'red', 'purple', 'orange'])
    ax.set_xlabel('Molécula')
    ax.set_ylabel('Intensidad Raman Promedio')
    ax.set_title('Comparación de Intensidad Raman entre Moléculas')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Añadir valores en las barras
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 2,
                f'{height}', ha='center', va='bottom')
    
    st.pyplot(fig)

def render_shift_estimation():
    st.markdown('<h2 class="section-header">Estimación de Desplazamiento Químico</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="highlight">
    <strong>Comparación de coordenadas ORCA</strong>
    </div>
    """, unsafe_allow_html=True)
    
    # Datos de ejemplo para la tabla de coordenadas
    coord_data = {
        "Molécula": ["ORCA"] * 10,
        "Átomo": ["C", "C", "C", "C", "C", "C", "C", "C", "N", "O"],
        "X (Å)": np.random.uniform(-2, 2, 10),
        "Y (Å)": np.random.uniform(-2, 2, 10),
        "Z (Å)": np.random.uniform(-2, 2, 10)
    }
    
    df = pd.DataFrame(coord_data)
    st.table(df)
    
    # Selector para NMR
    nmr_type = st.radio("Tipo de NMR:", ["¹H NMR", "¹³C NMR"])
    
    # Datos de ejemplo para NMR
    if nmr_type == "¹H NMR":
        shifts = np.random.uniform(0, 10, 5)
        atoms = ["H1", "H2", "H3", "H4", "H5"]
    else:
        shifts = np.random.uniform(0, 200, 5)
        atoms = ["C1", "C2", "C3", "C4", "C5"]
    
    # Crear gráfico de NMR
    fig, ax = plt.subplots(figsize=(10, 4))
    for i, (shift, atom) in enumerate(zip(shifts, atoms)):
        ax.plot([shift, shift], [0, 1], 'b-', linewidth=2)
        ax.text(shift, 1.05, atom, ha='center', va='bottom')
    
    ax.set_xlabel('Desplazamiento Químico (ppm)')
    ax.set_title(f'Espetro {nmr_type} Simulado')
    ax.set_yticks([])
    if nmr_type == "¹H NMR":
        ax.set_xlim(10, 0)  # NMR de hidrógeno va de 10 a 0 ppm
    else:
        ax.set_xlim(200, 0)  # NMR de carbono va de 200 a 0 ppm
    ax.grid(True, alpha=0.3, axis='x')
    
    st.pyplot(fig)

# Navegación principal
if section == "Molécula 3D":
    render_3d_molecule()
elif section == "Molécula 2D":
    render_2d_molecule()
elif section == "Conjunto de Moléculas":
    render_molecule_set()
elif section == "Contenedor de Moléculas":
    render_molecule_container()
elif section == "Espectro IR":
    render_ir_spectrum()
elif section == "Trabajo de Adhesión":
    render_adhesion_work()
elif section == "Molécula Teórica (RDF)":
    render_theoretical_molecule()
elif section == "Espectro Raman":
    render_raman_spectrum()
elif section == "Comparación de Moléculas":
    render_molecule_comparison()
elif section == "Raman vs Molécula":
    render_raman_vs_molecule()
elif section == "Estimación de Desplazamiento":
    render_shift_estimation()

# Pie de página
st.markdown("---")
st.caption("""
Sistema de Visualización Molecular - Desarrollado con Streamlit 🧪  
Basado en: Parra Figueredo, J. G., González, J., Perozo, E., & Iza, P. (2023).  
*Educación Química, 34*(1), 20-32.
""")

# Nota: Para las combinaciones necesarias en render_molecule_container
