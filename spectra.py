# spectra.py
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

# Importar nuestros parsers consolidados
from orca_parsers import (
    parse_orca_coordinates,
    parse_orca_frequencies,
    parse_orca_scf_energies,
    parse_orca_orbital_energies,
    parse_orca_xyz_file,
    get_molecule_outputs,
    check_orca_outputs_exist
)

# Importar nuestras clases personalizadas
from generate_inp import OrcaInputGenerator
from generar_out import OrcaOutputGenerator
from molecular_config_manager import MolecularConfigManager

# ========== FUNCIONES PARA LEER ARCHIVOS ORCA .out ==========

#----- Gráfico del Modelo 3D (DATOS REALES DE ORCA) -----
def dibujar_molecula_3d(molecule_name):
    """Visualización 3D de la molécula usando coordenadas optimizadas de ORCA"""
    
    # Verificar si existen archivos de ORCA
    outputs = get_molecule_outputs(molecule_name)
    
    if not outputs['opt'].exists():
        st.error(f"❌ No se encontró archivo de optimización: {outputs['opt']}")
        st.info("💡 Ejecuta primero el procesamiento con ORCA para generar los archivos necesarios.")
        return None
    
    # Leer coordenadas optimizadas
    elements, coords = parse_orca_coordinates(outputs['opt'])
    
    if len(elements) == 0:
        st.error("❌ No se pudieron leer las coordenadas del archivo ORCA")
        return None
    
    # Colores por elemento
    element_colors = {
        'H': '#FFFFFF',   # Blanco
        'C': '#404040',   # Gris oscuro  
        'N': '#3050F8',   # Azul
        'O': '#FF0D0D',   # Rojo
        'F': '#90E050',   # Verde claro
        'P': '#FF8000',   # Naranja
        'S': '#FFFF30',   # Amarillo
        'Cl': '#1FF01F',  # Verde
        'Br': '#A62929',  # Marrón rojizo
        'I': '#940094',   # Púrpura
    }
    
    # Radios por elemento (van der Waals)
    element_radii = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.64,
        'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39
    }
    
    # Crear figura 3D
    fig = go.Figure()
    
    # Función para crear esfera
    def create_sphere(center, radius, color):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = radius * np.cos(u) * np.sin(v) + center[0]
        y = radius * np.sin(u) * np.sin(v) + center[1] 
        z = radius * np.cos(v) + center[2]
        
        return go.Surface(
            x=x, y=y, z=z,
            colorscale=[[0, color], [1, color]],
            showscale=False,
            hoverinfo='skip',
            opacity=0.9
        )
    
    # Agregar átomos
    for i, (element, coord) in enumerate(zip(elements, coords)):
        color = element_colors.get(element, '#808080')  # Gris por defecto
        radius = element_radii.get(element, 0.5)
        
        sphere = create_sphere(coord, radius, color)
        fig.add_trace(sphere)
        
        # Agregar etiqueta del elemento
        fig.add_trace(go.Scatter3d(
            x=[coord[0]], y=[coord[1]], z=[coord[2]],
            mode='text',
            text=[f"{element}{i+1}"],
            textposition="middle center",
            textfont=dict(size=12, color='black'),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Agregar enlaces (distancia < 1.8 Å entre átomos no-H, < 1.2 Å para H)
    def create_cylinder(p1, p2, radius=0.1):
        # Vector del enlace
        vec = p2 - p1
        length = np.linalg.norm(vec)
        if length == 0:
            return None
            
        # Crear cilindro
        t = np.linspace(0, 1, 10)
        points = np.outer(t, vec) + p1
        
        return go.Scatter3d(
            x=points[:, 0], y=points[:, 1], z=points[:, 2],
            mode='lines',
            line=dict(width=8, color='gray'),
            showlegend=False,
            hoverinfo='skip'
        )
    
    # Detectar y dibujar enlaces
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            dist = np.linalg.norm(coords[i] - coords[j])
            
            # Criterios de enlace
            max_dist = 1.8
            if elements[i] == 'H' or elements[j] == 'H':
                max_dist = 1.2
                
            if dist < max_dist:
                bond = create_cylinder(coords[i], coords[j])
                if bond:
                    fig.add_trace(bond)
    
    # Configurar layout
    fig.update_layout(
        title=f"Molécula: {molecule_name} (Optimizada con ORCA)",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False), 
            zaxis=dict(visible=False),
            aspectmode='data',
            bgcolor='rgb(240, 240, 240)',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),
                center=dict(x=0, y=0, z=0)
            )
        ),
        margin=dict(l=0, r=0, t=50, b=0),
        showlegend=False,
        height=600
    )
    
    return fig

#----- Espectro IR (DATOS REALES DE ORCA) -----
def dibujar_espectro_ir(molecule_name):
    """Genera espectro IR usando datos reales de ORCA"""
    
    # Verificar archivos de ORCA
    outputs = get_molecule_outputs(molecule_name)
    
    if not outputs['ir-raman'].exists():
        st.error(f"❌ No se encontró archivo IR/Raman: {outputs['ir-raman']}")
        st.info("💡 Ejecuta primero el cálculo IR/Raman con ORCA.")
        return None
    
    # Leer datos de frecuencias
    spectra_data = parse_orca_frequencies(outputs['ir-raman'])
    ir_data = spectra_data['ir']
    
    if ir_data.empty:
        st.warning("⚠️ No se encontraron datos IR en el archivo de ORCA")
        return None
    
    # Crear el espectro
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if 'Frequency' in ir_data.columns and 'Intensity' in ir_data.columns:
        # Espectro de líneas
        frequencies = ir_data['Frequency'].values
        intensities = ir_data['Intensity'].values
        
        # Crear espectro gaussiano suavizado
        freq_range = np.linspace(max(frequencies) + 200, 400, 2000)  # Rango típico IR
        spectrum = np.zeros_like(freq_range)
        
        for freq, intensity in zip(frequencies, intensities):
            # Gaussiano centrado en cada frecuencia
            gaussian = intensity * np.exp(-0.5 * ((freq_range - freq) / 15) ** 2)
            spectrum += gaussian
        
        # Graficar espectro suavizado
        ax.plot(freq_range, spectrum, 'b-', linewidth=1.5, label=f'IR {molecule_name}')
        ax.fill_between(freq_range, spectrum, alpha=0.3)
        
        # Marcar picos principales
        for freq, intensity in zip(frequencies, intensities):
            if intensity > np.max(intensities) * 0.1:  # Solo picos significativos
                ax.axvline(x=freq, color='red', linestyle='--', alpha=0.7, linewidth=1)
                ax.text(freq, intensity * 1.1, f'{freq:.0f}', 
                       rotation=90, ha='center', va='bottom', fontsize=8)
    else:
        st.error("❌ Datos IR incompletos en el archivo ORCA")
        return None
    
    # Configurar gráfico
    ax.invert_xaxis()  # Típico en espectros IR
    ax.set_title(f"Espectro IR - {molecule_name} (calculado con ORCA)")
    ax.set_xlabel("Número de onda (cm⁻¹)")
    ax.set_ylabel("Intensidad (km/mol)")
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Mostrar información adicional
    st.subheader("📊 Información del espectro")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Modos vibracionales", len(ir_data))
    with col2:
        if 'Intensity' in ir_data.columns:
            st.metric("Intensidad máxima", f"{ir_data['Intensity'].max():.2f} km/mol")
    with col3:
        if 'Frequency' in ir_data.columns:
            st.metric("Frecuencia máxima", f"{ir_data['Frequency'].max():.0f} cm⁻¹")
    
    # Mostrar tabla de datos
    if st.checkbox("Ver tabla de frecuencias"):
        st.dataframe(ir_data, use_container_width=True)
    
    return fig

#----- Gráfico de Energías SCF (DATOS REALES DE ORCA) -----
def dibujar_energias_scf(molecule_name):
    """Genera gráficos de análisis energético usando datos reales de ORCA"""
    
    # Verificar archivos de ORCA
    outputs = get_molecule_outputs(molecule_name)
    
    if not outputs['opt'].exists():
        st.error(f"⚠️ No se encontró el archivo de optimización para {molecule_name}")
        st.info("💡 Ejecuta primero el procesamiento completo con ORCA para generar los archivos necesarios.")
        return None
    
    # Leer datos de energías SCF
    energias_eh = parse_orca_scf_energies(outputs['opt'])
    
    if not energias_eh:
        st.error(f"⚠️ No se pudieron extraer las energías SCF de {molecule_name}")
        return None
    
    # Crear figura con 2 gráficos (1x2)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(f'Análisis Energético - {molecule_name}', fontsize=16, fontweight='bold')
    
    # 1. Gráfico de energías electrónicas (Un electrón vs Dos electrones)
    electron_categorias = ['Nuclear\nRepulsion', 'One\nElectron', 'Two\nElectron', 'Electronic\nTotal', 'Total\nEnergy']
    electron_valores_eh = [
        energias_eh.get('Nuclear Repulsion', 0),
        energias_eh.get('One Electron Energy', 0),
        energias_eh.get('Two Electron Energy', 0),
        energias_eh.get('Electronic Energy', 0),
        energias_eh.get('Total Energy', 0)
    ]
    
    colors = ['#FF6B6B', '#8E44AD', '#E67E22', '#4ECDC4', '#2E86AB']
    bars = ax1.bar(electron_categorias, electron_valores_eh, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax1.set_ylabel('Energía (Eh)', fontweight='bold')
    ax1.set_title('Componentes de Energía Molecular (Hartree)', fontweight='bold')
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Añadir valores en las barras (Eh)
    for bar, valor in zip(bars, electron_valores_eh):
        if valor != 0:  # Solo mostrar si el valor existe
            height = bar.get_height()
            va_position = 'bottom' if height > 0 else 'top'
            y_offset = abs(height) * 0.02 if height > 0 else -abs(height) * 0.02
            ax1.text(bar.get_x() + bar.get_width()/2., height + y_offset,
                     f'{valor:.2f} Eh', ha='center', va=va_position, fontsize=8, fontweight='bold')
    
    # 2. Gráfico del teorema del virial
    virial_categorias = ['Energía Potencial', 'Energía Cinética', 'Energía Total']
    virial_valores_eh = [
        energias_eh.get('Potential Energy', 0),
        energias_eh.get('Kinetic Energy', 0),
        energias_eh.get('Total Energy', 0)
    ]
    
    bars_virial = ax2.bar(virial_categorias, virial_valores_eh, 
                         color=['#FF6B6B', '#4ECDC4', '#2E86AB'], alpha=0.8, 
                         edgecolor='black', linewidth=0.5)
    
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax2.set_ylabel('Energía (Eh)', fontweight='bold')
    ax2.set_title('Teorema del Virial', fontweight='bold')
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Añadir valores en las barras del virial
    for bar, valor in zip(bars_virial, virial_valores_eh):
        if valor != 0:  # Solo mostrar si el valor existe
            height = bar.get_height()
            va_position = 'bottom' if height > 0 else 'top'
            y_offset = abs(height) * 0.05 if height > 0 else -abs(height) * 0.05
            ax2.text(bar.get_x() + bar.get_width()/2., height + y_offset,
                     f'{valor:.2f} Eh', ha='center', va=va_position, fontsize=9, fontweight='bold')
    
    # Añadir ratio virial si está disponible
    if 'Virial Ratio' in energias_eh:
        virial_ratio = energias_eh['Virial Ratio']
        ax2.text(0.95, 0.95, f'Ratio Virial: {virial_ratio:.4f}\n(ideal ≈ 2.0)',
                 transform=ax2.transAxes, ha='right', va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=10)
    elif 'Potential Energy' in energias_eh and 'Kinetic Energy' in energias_eh and energias_eh['Kinetic Energy'] != 0:
        virial_ratio = abs(energias_eh['Potential Energy'] / energias_eh['Kinetic Energy'])
        ax2.text(0.95, 0.95, f'Ratio Virial: {virial_ratio:.4f}\n(ideal ≈ 2.0)',
                 transform=ax2.transAxes, ha='right', va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=10)
    
    plt.tight_layout()
    
    # Mostrar información adicional mejorada
    st.subheader("📊 Información energética detallada")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        if 'Total Energy' in energias_eh:
            st.metric("Energía Total", f"{energias_eh['Total Energy']:.6f} Eh")
    with col2:
        if 'One Electron Energy' in energias_eh:
            st.metric("Un Electrón", f"{energias_eh['One Electron Energy']:.6f} Eh")
    with col3:
        if 'Two Electron Energy' in energias_eh:
            st.metric("Dos Electrones", f"{energias_eh['Two Electron Energy']:.6f} Eh")
    with col4:
        if 'Virial Ratio' in energias_eh:
            st.metric("Ratio Virial", f"{energias_eh['Virial Ratio']:.6f}")
    
    # Información adicional sobre las energías electrónicas
    if all(key in energias_eh for key in ['One Electron Energy', 'Two Electron Energy', 'Electronic Energy']):
        st.subheader("⚡ Análisis de Energías Electrónicas")
        
        one_e = energias_eh['One Electron Energy']
        two_e = energias_eh['Two Electron Energy']
        total_e = energias_eh['Electronic Energy']
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.info(f"""
            **Energía de Un Electrón**
            - Valor: {one_e:.6f} Eh
            - Incluye: energía cinética de electrones y atracción núcleo-electrón
            - Siempre negativa (enlace estabilizante)
            """)
        
        with col2:
            st.success(f"""
            **Energía de Dos Electrones**
            - Valor: {two_e:.6f} Eh
            - Incluye: repulsión electrón-electrón
            - Siempre positiva (desestabilizante)
            """)
        
        with col3:
            balance = "Estabilizante" if total_e < 0 else "Desestabilizante"
            st.warning(f"""
            **Balance Electrónico**
            - Total: {total_e:.6f} Eh
            - Efecto neto: {balance}
            - Suma: {one_e:.3f} + {two_e:.3f} = {total_e:.3f}
            """)
    
    # Mostrar resumen numérico completo
    if st.checkbox("Ver resumen energético completo"):
        st.subheader("📋 Resumen Energético Completo")
        
        # Crear DataFrame para mostrar los datos
        data_rows = []
        for key, value in energias_eh.items():
            if not key.endswith('(eV)'):  # Solo mostrar valores en Hartree
                tipo = "Nuclear" if "Nuclear" in key else "Electrónica" if "Electronic" in key or "Electron" in key else "Mixta"
                ev_key = key + " (eV)"
                ev_value = energias_eh.get(ev_key, "N/A")
                data_rows.append({
                    'Componente': key,
                    'Energía (Eh)': f"{value:.6f}",
                    'Energía (eV)': f"{ev_value:.6f}" if ev_value != "N/A" else "N/A",
                    'Tipo': tipo
                })
        
        if data_rows:
            df_energias = pd.DataFrame(data_rows)
            st.dataframe(df_energias, use_container_width=True)
    
    return fig

#----- Gráfico de Energías Orbitales (DATOS REALES DE ORCA) -----
def dibujar_energias_orbitales(molecule_name):
    """Genera gráficos de análisis de energías orbitales usando datos reales de ORCA"""
    
    # Verificar archivos de ORCA
    outputs = get_molecule_outputs(molecule_name)
    
    if not outputs['opt'].exists():
        st.error(f"⚠️ No se encontró el archivo de optimización para {molecule_name}")
        return None
    
    # Leer datos de energías orbitales
    df = parse_orca_orbital_energies(outputs['opt'])
    
    if df.empty:
        st.error(f"⚠️ No se pudieron extraer las energías orbitales de {molecule_name}")
        return None
    
    # Clasificar orbitales
    df['Tipo'] = df['OCC'].apply(lambda x: 'Ocupado' if x > 0 else 'Virtual')
    df['HOMO_LUMO'] = 'Normal'
    
    # Identificar HOMO y LUMO
    ocupados = df[df['OCC'] > 0]
    virtuales = df[df['OCC'] == 0]
    
    homo_energy = lumo_energy = gap_energy = None
    
    if not ocupados.empty and not virtuales.empty:
        homo_idx = ocupados['E(Eh)'].idxmax()  # Mayor energía entre ocupados
        lumo_idx = virtuales['E(Eh)'].idxmin()  # Menor energía entre virtuales
        
        df.loc[homo_idx, 'HOMO_LUMO'] = 'HOMO'
        df.loc[lumo_idx, 'HOMO_LUMO'] = 'LUMO'
        
        homo_energy = df.loc[homo_idx, 'E(Eh)']
        lumo_energy = df.loc[lumo_idx, 'E(Eh)']
        gap_energy = lumo_energy - homo_energy
    
    # Crear figura con 3 subplots en filas diferentes (3x1)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15))
    fig.suptitle(f'Análisis de Energías Orbitales - {molecule_name} (Hartree)', fontsize=16, fontweight='bold')
    
    # 1. Diagrama de niveles de energía
    colors = {'Ocupado': 'blue', 'Virtual': 'red', 'HOMO': 'green', 'LUMO': 'orange'}
    
    for _, row in df.iterrows():
        color = colors[row['HOMO_LUMO']] if row['HOMO_LUMO'] != 'Normal' else colors[row['Tipo']]
        alpha = 1.0 if row['HOMO_LUMO'] != 'Normal' else 0.7
        
        ax1.hlines(y=row['NO'], xmin=df['E(Eh)'].min()-0.1, xmax=row['E(Eh)'], 
                  color=color, alpha=alpha, linewidth=3)
        ax1.plot(row['E(Eh)'], row['NO'], 'o', color=color, markersize=8, alpha=alpha)
        
        # Añadir etiqueta con valor energético para HOMO y LUMO
        if row['HOMO_LUMO'] != 'Normal':
            ax1.text(row['E(Eh)'] + 0.02, row['NO'], f'{row["E(Eh)"]:.3f}', 
                    va='center', fontweight='bold', fontsize=9)
    
    ax1.set_xlabel('Energía (Eh)', fontweight='bold')
    ax1.set_ylabel('Número Orbital', fontweight='bold')
    ax1.set_title('Diagrama de Niveles de Energía', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(df['E(Eh)'].min()-0.1, df['E(Eh)'].max()+0.1)
    ax1.set_ylim(-0.5, len(df)+0.5)
    
    # Añadir línea en energía cero
    ax1.axvline(x=0, color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    # 2. Diagrama de dispersión Ocupación vs Energía
    scatter = ax2.scatter(df['E(Eh)'], df['OCC'], 
                         c=df['E(Eh)'], cmap='viridis', 
                         s=80, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='E=0 Eh')
    ax2.set_xlabel('Energía (Eh)', fontweight='bold')
    ax2.set_ylabel('Ocupación', fontweight='bold')
    ax2.set_title('Ocupación vs Energía', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax2, label='Energía (Eh)')
    
    # 3. Gráfico de barras horizontal de energías orbitales
    colors_bars = []
    for _, row in df.iterrows():
        if row['HOMO_LUMO'] == 'HOMO':
            colors_bars.append('green')
        elif row['HOMO_LUMO'] == 'LUMO':
            colors_bars.append('orange')
        elif row['Tipo'] == 'Ocupado':
            colors_bars.append('blue')
        else:
            colors_bars.append('red')
    
    # Usar barh para barras horizontales
    bars = ax3.barh(range(len(df)), df['E(Eh)'], color=colors_bars, alpha=0.8, edgecolor='black')
    ax3.axvline(x=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax3.set_ylabel('Número Orbital', fontweight='bold')
    ax3.set_xlabel('Energía (Eh)', fontweight='bold')
    ax3.set_title('Energías Orbitales Individuales (Horizontal)', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='x')
    
    # Añadir etiquetas a las barras horizontales solo para valores significativos
    for i, (bar, valor) in enumerate(zip(bars, df['E(Eh)'])):
        if abs(valor) > 0.1:  # Solo etiquetar valores significativos
            width = bar.get_width()
            ha = 'left' if width > 0 else 'right'
            x_offset = 0.02 if width > 0 else -0.02
            ax3.text(width + x_offset, bar.get_y() + bar.get_height()/2.,
                    f'{valor:.3f}', ha=ha, va='center', fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    
    # Mostrar información adicional mejorada
    st.subheader("🔬 Información de Energías Orbitales")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Orbitales", len(df))
    with col2:
        st.metric("Ocupados", len(ocupados))
    with col3:
        st.metric("Virtuales", len(virtuales))
    with col4:
        if gap_energy is not None:
            st.metric("Gap HOMO-LUMO", f"{gap_energy:.6f} Eh")
    
    # Información sobre HOMO y LUMO
    if homo_energy is not None and lumo_energy is not None:
        st.subheader("⚡ Análisis HOMO-LUMO")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.info(f"""
            **HOMO (Orbital Molecular Ocupado más Alto)**
            - Energía: {homo_energy:.6f} Eh
            - Energía: {homo_energy*27.2114:.3f} eV
            - Capacidad donadora de electrones
            """)
        
        with col2:
            st.success(f"""
            **LUMO (Orbital Molecular No Ocupado más Bajo)**
            - Energía: {lumo_energy:.6f} Eh
            - Energía: {lumo_energy*27.2114:.3f} eV
            - Capacidad aceptora de electrones
            """)
        
        with col3:
            st.warning(f"""
            **Gap HOMO-LUMO**
            - Diferencia: {gap_energy:.6f} Eh
            - Diferencia: {gap_energy*27.2114:.3f} eV
            - Indicador de reactividad química
            """)
    
    # Estadísticas detalladas
    if st.checkbox("Ver estadísticas orbitales detalladas"):
        st.subheader("📊 Estadísticas de Energías Orbitales")
        
        # Crear tabla de estadísticas
        stats_data = []
        
        if not ocupados.empty:
            stats_data.extend([
                {"Categoría": "Ocupados - Energía Mínima", "Valor": f"{ocupados['E(Eh)'].min():.6f} Eh", "eV": f"{ocupados['E(Eh)'].min()*27.2114:.3f} eV"},
                {"Categoría": "Ocupados - Energía Máxima", "Valor": f"{ocupados['E(Eh)'].max():.6f} Eh", "eV": f"{ocupados['E(Eh)'].max()*27.2114:.3f} eV"},
                {"Categoría": "Ocupados - Promedio", "Valor": f"{ocupados['E(Eh)'].mean():.6f} Eh", "eV": f"{ocupados['E(Eh)'].mean()*27.2114:.3f} eV"},
            ])
        
        if not virtuales.empty:
            stats_data.extend([
                {"Categoría": "Virtuales - Energía Mínima", "Valor": f"{virtuales['E(Eh)'].min():.6f} Eh", "eV": f"{virtuales['E(Eh)'].min()*27.2114:.3f} eV"},
                {"Categoría": "Virtuales - Energía Máxima", "Valor": f"{virtuales['E(Eh)'].max():.6f} Eh", "eV": f"{virtuales['E(Eh)'].max()*27.2114:.3f} eV"},
                {"Categoría": "Virtuales - Promedio", "Valor": f"{virtuales['E(Eh)'].mean():.6f} Eh", "eV": f"{virtuales['E(Eh)'].mean()*27.2114:.3f} eV"},
            ])
        
        if stats_data:
            df_stats = pd.DataFrame(stats_data)
            st.dataframe(df_stats, use_container_width=True)
    
    # Mostrar tabla de orbitales críticos
    if st.checkbox("Ver tabla de orbitales críticos"):
        st.subheader("🎯 Orbitales Críticos (HOMO y LUMO)")
        critical_orbitals = df[df['HOMO_LUMO'].isin(['HOMO', 'LUMO'])]
        if not critical_orbitals.empty:
            critical_display = critical_orbitals[['NO', 'OCC', 'E(Eh)', 'HOMO_LUMO']].copy()
            critical_display['E(eV)'] = critical_display['E(Eh)'] * 27.2114
            st.dataframe(critical_display, use_container_width=True)
    
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
        vecs[:, 2] += lift  # elevar un poco en Z
        return vecs
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

#----- Conjunto de moléculas dinámico (usando datos reales de ORCA) -----
def dibujar_conjunto_molecula(molecule_name):
    """Visualización de conjunto de moléculas usando coordenadas optimizadas de ORCA"""
    
    # Verificar si existen archivos de ORCA
    outputs = get_molecule_outputs(molecule_name)
    
    if not outputs['opt'].exists():
        st.error(f"❌ No se encontró archivo de optimización: {outputs['opt']}")
        return None
    
    # Leer coordenadas optimizadas
    elements, coords = parse_orca_coordinates(outputs['opt'])
    
    if len(elements) == 0:
        st.error("❌ No se pudieron leer las coordenadas del archivo ORCA")
        return None
    
    # Colores por elemento
    element_colors = {
        'H': '#FFFFFF',   # Blanco
        'C': '#404040',   # Gris oscuro  
        'N': '#3050F8',   # Azul
        'O': '#FF0D0D',   # Rojo
        'F': '#90E050',   # Verde claro
        'P': '#FF8000',   # Naranja
        'S': '#FFFF30',   # Amarillo
        'Cl': '#1FF01F',  # Verde
        'Br': '#A62929',  # Marrón rojizo
        'I': '#940094',   # Púrpura
    }
    
    # Radios por elemento
    element_radii = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.64,
        'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39
    }
    
    # Helpers geométricos
    def create_sphere(center, radius, color):
        u, v = np.mgrid[0:2*np.pi:12j, 0:np.pi:6j]
        x = radius * np.cos(u) * np.sin(v) + center[0]
        y = radius * np.sin(u) * np.sin(v) + center[1]
        z = radius * np.cos(v) + center[2]
        return go.Surface(
            x=x, y=y, z=z,
            colorscale=[[0, color], [1, color]],
            showscale=False, hoverinfo='skip', opacity=0.9
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

    # Crear figura
    fig = go.Figure()
    
    # Posiciones para el conjunto de moléculas (arreglo lineal conectado)
    spacing = max(np.max(coords, axis=0) - np.min(coords, axis=0)) * 2.0  # Espaciado dinámico
    
    # Crear cadena de moléculas con patrón zig-zag como en NH3
    num_molecules = 7  # Número de moléculas en la cadena
    grid_positions = []
    
    for i in range(num_molecules):
        # Patrón zig-zag más pronunciado alternando en X y Z
        x_offset = 0.8 if i % 2 == 0 else -0.8  # Aumentado de 0.3 a 0.8
        y_offset = (num_molecules - 1 - i) * spacing  # De arriba hacia abajo
        z_offset = 0.4 if i % 2 == 0 else -0.4  # Añadido zig-zag en Z también
        
        offset = np.array([x_offset, y_offset, z_offset])
        grid_positions.append(offset)
    
    # Calcular centros de masa de cada molécula para las conexiones
    molecule_centers = []
    for grid_pos in grid_positions:
        # Centro de masa de la molécula
        center_of_mass = np.mean(coords, axis=0) + grid_pos
        molecule_centers.append(center_of_mass)
    
    # Dibujar cada molécula en las posiciones del grid
    for grid_pos in grid_positions:
        # Agregar átomos
        for element, coord in zip(elements, coords):
            new_pos = coord + grid_pos
            color = element_colors.get(element, '#808080')
            radius = element_radii.get(element, 0.5) * 0.8  # Reducir un poco para el conjunto
            
            sphere = create_sphere(new_pos, radius, color)
            fig.add_trace(sphere)
        
        # Agregar enlaces dentro de cada molécula
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                
                # Criterios de enlace
                max_dist = 1.8
                if elements[i] == 'H' or elements[j] == 'H':
                    max_dist = 1.2
                    
                if dist < max_dist:
                    pos1 = coords[i] + grid_pos
                    pos2 = coords[j] + grid_pos
                    bond = create_cylinder(pos1, pos2, 0.05, '#808080')
                    if bond:
                        fig.add_trace(bond)
    
    # Agregar conexiones entre moléculas (enlaces intermoleculares con zig-zag)
    connection_color = '#FFFFFF'  # Color blanco para las conexiones
    for i in range(len(molecule_centers) - 1):
        center1 = molecule_centers[i]
        center2 = molecule_centers[i + 1]
        
        # Crear conexión zig-zag más suave con punto intermedio
        # Punto intermedio para crear curva zig-zag
        mid_point = (center1 + center2) / 2
        
        # Añadir desplazamiento perpendicular para hacer la curva más visible
        direction = center2 - center1
        perpendicular = np.array([-direction[2], 0, direction[0]])  # Vector perpendicular
        if np.linalg.norm(perpendicular) > 0:
            perpendicular = perpendicular / np.linalg.norm(perpendicular)
            # Alternar la dirección de la curva
            curve_offset = 0.3 * (1 if i % 2 == 0 else -1)
            mid_point += perpendicular * curve_offset
        
        # Crear dos segmentos: center1 -> mid_point -> center2
        connection1 = create_cylinder(center1, mid_point, 0.08, connection_color)
        connection2 = create_cylinder(mid_point, center2, 0.08, connection_color)
        
        if connection1:
            fig.add_trace(connection1)
        if connection2:
            fig.add_trace(connection2)
        
        # Agregar una pequeña esfera en el punto medio para suavizar la conexión
        mid_sphere = create_sphere(mid_point, 0.06, connection_color)
        fig.add_trace(mid_sphere)
    
    # Configurar layout
    fig.update_layout(
        title=f"Cadena conectada de moléculas: {molecule_name} (7 moléculas enlazadas)",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False), 
            zaxis=dict(visible=False),
            aspectmode='data',
            bgcolor='rgb(14,17,23)',  # Fondo oscuro como en NH3
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),
                center=dict(x=0, y=0, z=0)
            )
        ),
        margin=dict(l=0, r=0, t=50, b=0),
        showlegend=False,
        height=600
    )
    
    return fig

#----- Gráfico del Modelo 3D reactivo (usando configuración automática) -----
def dibujar_molecula_3d(molecule_name):
    """Visualización 3D de la molécula usando datos de configuración automática"""
    
    try:
        config_manager = MolecularConfigManager()
        elements, coords, description = config_manager.read_xyz_file(molecule_name)
        
        if len(elements) == 0:
            st.error("❌ No se pudieron leer los datos de la molécula")
            return None
            
    except Exception as e:
        st.error(f"❌ Error leyendo datos de la molécula: {e}")
        return None
    
    # Colores por elemento
    element_colors = {
        'H': '#FFFFFF',   # Blanco
        'C': '#404040',   # Gris oscuro  
        'N': '#3050F8',   # Azul
        'O': '#FF0D0D',   # Rojo
        'F': '#90E050',   # Verde claro
        'P': '#FF8000',   # Naranja
        'S': '#FFFF30',   # Amarillo
        'Cl': '#1FF01F',  # Verde
        'Br': '#A62929',  # Marrón rojizo
        'I': '#940094',   # Púrpura
    }
    
    # Radios por elemento
    element_radii = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.64,
        'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39
    }
    
    # Helpers geométricos
    def create_sphere(center, radius, color):
        u, v = np.mgrid[0:2*np.pi:12j, 0:np.pi:6j]
        x = radius * np.cos(u) * np.sin(v) + center[0]
        y = radius * np.sin(u) * np.sin(v) + center[1]
        z = radius * np.cos(v) + center[2]
        return go.Surface(
            x=x, y=y, z=z,
            colorscale=[[0, color], [1, color]],
            showscale=False, hoverinfo='skip', opacity=0.9
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

    # Crear figura
    fig = go.Figure()
    
    # Agregar átomos
    for element, coord in zip(elements, coords):
        color = element_colors.get(element, '#808080')
        radius = element_radii.get(element, 0.5)
        
        sphere = create_sphere(coord, radius, color)
        fig.add_trace(sphere)
    
    # Agregar enlaces 
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            dist = np.linalg.norm(coords[i] - coords[j])
            
            # Criterios de enlace
            max_dist = 1.8
            if elements[i] == 'H' or elements[j] == 'H':
                max_dist = 1.2
                
            if dist < max_dist:
                bond = create_cylinder(coords[i], coords[j], 0.08, '#808080')
                if bond:
                    fig.add_trace(bond)
    
    # Configurar layout
    fig.update_layout(
        title=f"Molécula 3D: {molecule_name}",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False), 
            zaxis=dict(visible=False),
            aspectmode='data',
            bgcolor='rgb(240,240,240)',  # Fondo claro
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),
                center=dict(x=0, y=0, z=0)
            )
        ),
        margin=dict(l=0, r=0, t=50, b=0),
        showlegend=False,
        height=600
    )
    
    # Mostrar información adicional
    st.subheader("📊 Información de la molécula")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Número de átomos", len(elements))
    with col2:
        elementos_unicos = list(set(elements))
        st.metric("Tipos de elementos", len(elementos_unicos))
    with col3:
        # Calcular número de enlaces
        num_enlaces = 0
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                max_dist = 1.8
                if elements[i] == 'H' or elements[j] == 'H':
                    max_dist = 1.2
                if dist < max_dist:
                    num_enlaces += 1
        st.metric("Enlaces detectados", num_enlaces)
    
    if description and description.strip():
        st.info(f"**Descripción:** {description}")
    
    return fig

#----- Gráfico del Modelo 2D (DATOS REALES DE ORCA) -----
def dibujar_molecula_2d(molecule_name):
    """Visualización 2D de la molécula usando coordenadas optimizadas de ORCA"""
    
    # Verificar si existen archivos de ORCA
    outputs = get_molecule_outputs(molecule_name)
    
    if not outputs['opt'].exists():
        st.error(f"❌ No se encontró archivo de optimización: {outputs['opt']}")
        st.info("💡 Ejecuta primero el procesamiento con ORCA para generar los archivos necesarios.")
        return None
    
    # Leer coordenadas optimizadas
    elements, coords = parse_orca_coordinates(outputs['opt'])
    
    if len(elements) == 0:
        st.error("❌ No se pudieron leer las coordenadas del archivo ORCA")
        return None
    
    # Colores por elemento (más brillantes para 2D)
    element_colors = {
        'H': '#E0E0E0',   # Gris claro
        'C': '#202020',   # Negro
        'N': '#4070FF',   # Azul brillante
        'O': '#FF4040',   # Rojo brillante
        'F': '#90FF50',   # Verde brillante
        'P': '#FFA040',   # Naranja
        'S': '#FFFF40',   # Amarillo
        'Cl': '#40FF40',  # Verde
        'Br': '#C04040',  # Marrón rojizo
        'I': '#B040B0',   # Púrpura
    }
    
    # Radios por elemento (escalados para 2D)
    element_radii = {
        'H': 0.4, 'C': 0.8, 'N': 0.75, 'O': 0.7, 'F': 0.65,
        'P': 1.1, 'S': 1.0, 'Cl': 1.05, 'Br': 1.25, 'I': 1.4
    }
    
    # Proyectar a 2D (usar coordenadas X, Y)
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    
    # Crear figura
    fig, ax = plt.subplots(figsize=(10, 8), dpi=120)
    fig.patch.set_facecolor('#FAFAFA')
    ax.set_facecolor('#FAFAFA')
    
    # Función para dibujar enlaces
    def draw_bond(ax, p1, p2, bond_width=6):
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 
                color='#666666', linewidth=bond_width, 
                solid_capstyle='round', zorder=1)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 
                color='#CCCCCC', linewidth=bond_width-2, 
                solid_capstyle='round', zorder=1.1)
    
    # Dibujar enlaces (distancia < 1.8 Å)
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            dist = np.linalg.norm(coords[i] - coords[j])
            
            max_dist = 1.8
            if elements[i] == 'H' or elements[j] == 'H':
                max_dist = 1.2
                
            if dist < max_dist:
                draw_bond(ax, coords[i][:2], coords[j][:2])
    
    # Función para dibujar átomos con efecto 3D
    def draw_glossy_atom(ax, center, radius, face_color, edge_color):
        # Sombra
        shadow = plt.Circle((center[0]+0.15, center[1]-0.15), radius*1.05,
                           facecolor='black', edgecolor='none', 
                           alpha=0.2, zorder=0.5)
        ax.add_patch(shadow)
        
        # Átomo principal
        atom = plt.Circle(center, radius, 
                         facecolor=face_color, edgecolor=edge_color,
                         linewidth=3, zorder=2)
        ax.add_patch(atom)
        
        # Brillos
        highlight1 = plt.Circle((center[0]-radius*0.3, center[1]+radius*0.3), 
                               radius*0.4, facecolor='white', 
                               edgecolor='none', alpha=0.6, zorder=3)
        ax.add_patch(highlight1)
        
        highlight2 = plt.Circle((center[0]-radius*0.15, center[1]+radius*0.15), 
                               radius*0.2, facecolor='white', 
                               edgecolor='none', alpha=0.8, zorder=3)
        ax.add_patch(highlight2)
    
    # Dibujar átomos
    for i, (element, coord) in enumerate(zip(elements, coords)):
        face_color = element_colors.get(element, '#808080')
        edge_color = face_color.replace('FF', 'CC').replace('40', '20')  # Color más oscuro
        radius = element_radii.get(element, 0.6)
        
        draw_glossy_atom(ax, coord[:2], radius, face_color, edge_color)
        
        # Etiqueta del átomo
        ax.text(coord[0], coord[1], f"{element}{i+1}",
               ha='center', va='center', fontsize=12, fontweight='bold',
               color='white' if element in ['C', 'N'] else 'black',
               zorder=4)
    
    # Configurar ejes
    margin = 2.0
    ax.set_xlim(x_coords.min() - margin, x_coords.max() + margin)
    ax.set_ylim(y_coords.min() - margin, y_coords.max() + margin)
    ax.set_aspect('equal', adjustable='box')
    
    # Título y etiquetas
    ax.set_title(f"Molécula: {molecule_name} (vista 2D)", fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel("X (Angstrom)", fontsize=12)
    ax.set_ylabel("Y (Angstrom)", fontsize=12)
    
    # Grid sutil
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Información molecular en el gráfico
    info_text = f"Átomos: {len(elements)}\nElementos: {', '.join(set(elements))}"
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
           verticalalignment='top', bbox=dict(boxstyle='round', 
           facecolor='white', alpha=0.8), fontsize=10)
    
    return fig
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

    y_vdw_contacto = np.array([2, 0.85, 0.60, 0.08, 0.05], dtype=float)
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

# ---------- NUEVA FUNCIÓN PARA SELECTOR DE MOLÉCULAS ----------
def seleccionar_molecula():
    """Función para seleccionar una molécula de la carpeta moleculas_xyz."""
    moleculas_dir = Path("moleculas_xyz")
    
    if not moleculas_dir.exists():
        st.error("La carpeta 'moleculas_xyz' no existe")
        return None
        
    archivos_xyz = list(moleculas_dir.glob("*.xyz"))
    
    if not archivos_xyz:
        st.error("No se encontraron archivos .xyz en la carpeta 'moleculas_xyz'")
        return None
        
    nombres_moleculas = [f.stem for f in archivos_xyz]
    
    # Inicializar el gestor de configuración
    try:
        config_manager = MolecularConfigManager()
        
        # Detectar molécula actualmente configurada
        current_molecule = config_manager.get_current_molecule_from_paso(2)
        
        # Encontrar el índice de la molécula actual
        default_index = 0
        if current_molecule and current_molecule in nombres_moleculas:
            default_index = nombres_moleculas.index(current_molecule)
            st.sidebar.info(f"🎯 Molécula detectada: {current_molecule}")
        
    except Exception as e:
        st.sidebar.warning(f"⚠️ Error detectando molécula actual: {e}")
        default_index = 0
    
    # Selector de molécula
    molecula_seleccionada = st.selectbox(
        "Selecciona una molécula para analizar:",
        nombres_moleculas,
        index=default_index,
        key="molecule_selector"
    )
    
    # Sistema reactivo: configurar automáticamente cuando cambia la selección
    if molecula_seleccionada:
        try:
            # Verificar si necesita reconfiguración
            config_manager = MolecularConfigManager()
            current_configured = config_manager.get_current_molecule_from_paso(2)
            
            if current_configured != molecula_seleccionada:
                st.info(f"🔄 Configurando automáticamente: {molecula_seleccionada}")
                
                # Auto-configurar la nueva molécula
                with st.spinner(f"Actualizando archivos de configuración para {molecula_seleccionada}..."):
                    results = config_manager.auto_configure_molecule(molecula_seleccionada)
                
                # Mostrar resultados de la configuración
                successful_updates = sum(results.values())
                total_updates = len(results)
                
                if successful_updates == total_updates:
                    st.success(f"✅ {molecula_seleccionada} configurada correctamente")
                else:
                    st.warning(f"⚠️ Configuración parcial: {successful_updates}/{total_updates} archivos actualizados")
                
                # Mostrar detalles de cada paso
                with st.expander("Ver detalles de configuración"):
                    for paso, success in results.items():
                        status = "✅" if success else "❌"
                        st.write(f"{status} {paso}: {'Exitoso' if success else 'Falló'}")
                        
        except Exception as e:
            st.error(f"❌ Error en configuración automática: {e}")
    
    return molecula_seleccionada

def procesar_molecula_completa(nombre_molecula: str):
    """Función que ejecuta todo el pipeline: generar inputs -> ejecutar ORCA -> mostrar resultados."""
    
    if not nombre_molecula:
        st.warning("No se ha seleccionado ninguna molécula")
        return
        
    molecula_path = Path("moleculas_xyz") / f"{nombre_molecula}.xyz"
    
    if not molecula_path.exists():
        st.error(f"No se encontró el archivo: {molecula_path}")
        return
        
    # Mostrar información de la molécula
    st.info(f"Procesando molécula: {nombre_molecula}")
    
    # Crear columnas para organizar la información
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("🔄 Generando archivos de entrada ORCA")
        
        # Generar archivos .inp
        generator = OrcaInputGenerator(str(molecula_path))
        generator.generate_all()
        
        st.success("✅ Archivos .inp generados correctamente")
        
    with col2:
        st.subheader("🖥️ Ejecutando cálculos ORCA")
        
        # Ejecutar ORCA (si está disponible)
        output_generator = OrcaOutputGenerator()
        
        if output_generator.orca_path:
            if st.button("Ejecutar cálculos ORCA"):
                with st.spinner("Ejecutando cálculos... Esto puede tomar varios minutos"):
                    results = output_generator.process_molecule(nombre_molecula)
                    
                if results["success"]:
                    st.success("✅ Cálculos ORCA completados")
                    
                    # Mostrar resultados de cada cálculo
                    for calc_type, calc_result in results["calculations"].items():
                        if calc_result["success"]:
                            tiempo = calc_result["time"]
                            st.write(f"- {calc_type.upper()}: Completado en {tiempo:.1f}s")
                        else:
                            st.warning(f"- {calc_type.upper()}: Falló")
                else:
                    st.error("❌ Algunos cálculos fallaron")
        else:
            st.warning("⚠️ ORCA no está disponible. Los cálculos no se pueden ejecutar.")
            st.info("Los archivos .inp se han generado y pueden ser ejecutados manualmente.")

def mostrar_estado_configuracion():
    """Mostrar el estado actual de configuración de la aplicación."""
    st.subheader("⚙️ Estado de Configuración")
    
    try:
        config_manager = MolecularConfigManager()
        
        # Detectar molécula configurada
        current_molecule = config_manager.get_current_molecule_from_paso(2)
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if current_molecule:
                st.metric("Molécula Actual", current_molecule)
            else:
                st.metric("Molécula Actual", "No detectada")
        
        with col2:
            molecules = config_manager.get_available_molecules()
            st.metric("Moléculas Disponibles", len(molecules))
        
        with col3:
            # Verificar archivos de pasos
            pasos_existentes = 0
            for paso_num, paso_path in config_manager.paso_files.items():
                if paso_path.exists():
                    pasos_existentes += 1
            st.metric("Archivos de Paso", f"{pasos_existentes}/3")
        
        # Tabla de archivos de pasos
        st.subheader("📁 Estado de Archivos")
        
        data = []
        for paso_num, paso_path in config_manager.paso_files.items():
            exists = paso_path.exists()
            if exists:
                # Obtener tamaño del archivo
                size = paso_path.stat().st_size
                modified = paso_path.stat().st_mtime
                import time
                mod_time = time.ctime(modified)
                status = "✅ Existe"
                info = f"{size} bytes, modificado: {mod_time}"
            else:
                status = "❌ No existe"
                info = "Archivo no encontrado"
            
            data.append({
                "Paso": f"paso_{paso_num}.txt",
                "Estado": status,
                "Información": info
            })
        
        df_estado = pd.DataFrame(data)
        st.dataframe(df_estado, use_container_width=True)
        
        # Botón para reconfigurar todo
        st.subheader("🔄 Acciones")
        
        if current_molecule:
            if st.button(f"Reconfigurar {current_molecule}", type="primary"):
                with st.spinner("Reconfigurando..."):
                    results = config_manager.auto_configure_molecule(current_molecule)
                
                successful_updates = sum(results.values())
                total_updates = len(results)
                
                if successful_updates == total_updates:
                    st.success("✅ Reconfiguración completada")
                    st.rerun()
                else:
                    st.warning(f"⚠️ Reconfiguración parcial: {successful_updates}/{total_updates}")
        
        # Limpiar configuración
        if st.button("🗑️ Limpiar configuración", type="secondary"):
            if st.checkbox("Confirmar limpieza"):
                for paso_path in config_manager.paso_files.values():
                    try:
                        if paso_path.exists():
                            paso_path.unlink()
                        st.success(f"Eliminado: {paso_path.name}")
                    except Exception as e:
                        st.error(f"Error eliminando {paso_path.name}: {e}")
                
                st.rerun()
        
    except Exception as e:
        st.error(f"❌ Error obteniendo estado de configuración: {e}")

def mostrar_informacion_molecula(nombre_molecula: str):
    """Muestra información básica sobre la molécula seleccionada."""
    if not nombre_molecula:
        st.warning("No hay molécula seleccionada")
        return
        
    try:
        config_manager = MolecularConfigManager()
        elements, coordinates, description = config_manager.read_xyz_file(nombre_molecula)
        
        # Mostrar información
        st.subheader(f"📋 Información de {nombre_molecula}")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Número de átomos", len(elements))
            
        with col2:
            elementos_unicos = list(set(elements))
            st.metric("Tipos de elementos", len(elementos_unicos))
            
        with col3:
            formula = ""
            from collections import Counter
            contador = Counter(elements)
            for elemento in sorted(contador.keys()):
                count = contador[elemento]
                if count > 1:
                    formula += f"{elemento}₂" if count == 2 else f"{elemento}{count}"
                else:
                    formula += elemento
            st.metric("Fórmula", formula)
        
        # Descripción
        if description and description.strip():
            st.write(f"**Descripción:** {description}")
            
        # Tabla de coordenadas
        if len(elements) > 0:
            df_coords = pd.DataFrame({
                "Elemento": elements,
                "X (Å)": coordinates[:, 0],
                "Y (Å)": coordinates[:, 1], 
                "Z (Å)": coordinates[:, 2]
            })
            st.subheader("🧮 Coordenadas atómicas")
            st.dataframe(df_coords, use_container_width=True)
            
        # Información adicional
        st.subheader("📊 Análisis Geométrico")
        
        # Centro de masa
        center_of_mass = np.mean(coordinates, axis=0)
        st.write(f"**Centro de masa:** ({center_of_mass[0]:.3f}, {center_of_mass[1]:.3f}, {center_of_mass[2]:.3f}) Å")
        
        # Dimensiones de la caja
        min_coords = np.min(coordinates, axis=0)
        max_coords = np.max(coordinates, axis=0)
        dimensions = max_coords - min_coords
        st.write(f"**Dimensiones:** {dimensions[0]:.3f} × {dimensions[1]:.3f} × {dimensions[2]:.3f} Å")
        
        # Radio de giro
        if len(coordinates) > 1:
            distances_from_center = np.linalg.norm(coordinates - center_of_mass, axis=1)
            radius_of_gyration = np.sqrt(np.mean(distances_from_center**2))
            st.write(f"**Radio de giro:** {radius_of_gyration:.3f} Å")
        
    except Exception as e:
        st.error(f"❌ Error cargando información de la molécula: {e}")
