# main.py
import streamlit as st
from pathlib import Path
import os
import pandas as pd

from rdf_analysis import analyze_molecule_rdf
# Importar funciones de spectra
from spectra import (
    seleccionar_molecula,
    procesar_molecula_completa,
    mostrar_estado_configuracion,
    mostrar_informacion_molecula,
    dibujar_molecula_3d,
    dibujar_molecula_2d,
    dibujar_conjunto_molecula,
    contenedor_molecular,
    dibujar_espectro_ir,
    dibujar_energias_scf,
    dibujar_energias_orbitales,
    dibujar_analisis_poblacion,
    plot_ir_spectrum,
    plot_ir_spectrum_comparison,
    graficar_trabajo_adhesion,
    mostrar_rdf,
    mostrar_ir_raman,
    comparar_moleculas_orca_vs_nh3,
    comparar_ir_raman_vs_nh3,
    comparar_rmn_s4_vs_nh3,
    get_molecule_outputs
)

# Importar el gestor de configuración
from molecular_config_manager import MolecularConfigManager
from generar_out import parse_vibrational_frequencies

# ---------- MAIN MODIFICADO ----------
def main():
    st.set_page_config(
        page_title="Analizador de Moléculas XYZ",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 Analizador de Moléculas XYZ con ORCA")
    st.markdown("---")
    
    # Sidebar para navegación
    st.sidebar.title("🗂️ Navegación")
    
    # Primer paso: Selección de molécula
    st.sidebar.subheader("1. Selección de Molécula")
    molecula_seleccionada = seleccionar_molecula()
    
    if molecula_seleccionada:
        st.sidebar.success(f"Molécula seleccionada: {molecula_seleccionada}")
        
        # Opciones de análisis
        st.sidebar.subheader("2. Tipo de Análisis")
        option = st.sidebar.radio("Selecciona una sección:", [
            "⚙️ Estado de configuración",
            "📋 Información de la molécula",
            "🔄 Procesar con ORCA",
            "🧪 Molécula 3D",
            "📊 Molécula 2D", 
            "🔗 Conjunto de moléculas",
            "📦 Contenedor de moléculas",
            "📈 Espectro IR",
            "⚡ Energías SCF",
            "🔬 Energías Orbitales",
            "🧬 Análisis de Población (Mulliken/Löwdin)",
            "🔬 Trabajo de adhesión",
            "📊 Función de Distribución Radial",
            "📉 Espectro IR Teórico",
            "⚛️ Molécula teórica (RDF)",
            "🌈 Espectros Raman",
            "🔍 Comparación con NH₃",
            "🧬 IR/Raman vs NH₃",
            "📉 Desplazamientos químicos"
        ])
        
        # Contenido principal
        if option == "⚙️ Estado de configuración":
            st.header("⚙️ Estado de Configuración del Sistema")
            mostrar_estado_configuracion()
            
        elif option == "📋 Información de la molécula":
            mostrar_informacion_molecula(molecula_seleccionada)
            
        elif option == "🔄 Procesar con ORCA":
            st.header("🔄 Procesamiento completo con ORCA")
            st.write("Este proceso generará los archivos de entrada de ORCA y ejecutará los cálculos.")
            
            with st.expander("ℹ️ Información del proceso", expanded=True):
                st.write("""
                **El procesamiento incluye:**
                1. 📝 Generación de archivos .inp (optimización, IR/Raman, RMN)
                2. ⚙️ Ejecución de cálculos ORCA (si está disponible)
                3. 📊 Análisis de resultados
                
                **Nota:** Se requiere tener ORCA instalado para ejecutar los cálculos.
                """)
            
            procesar_molecula_completa(molecula_seleccionada)
            
        elif option == "🧪 Molécula 3D":
            st.header(f"🧪 Visualización 3D - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
                return
            
            fig = dibujar_molecula_3d(molecula_seleccionada)
            if fig:
                st.plotly_chart(fig, width='stretch')
            
        elif option == "📊 Molécula 2D":
            st.header(f"📊 Visualización 2D - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
                return
            
            fig = dibujar_molecula_2d(molecula_seleccionada)
            if fig:
                st.pyplot(fig)
            
        elif option == "🔗 Conjunto de moléculas":
            st.header(f"🔗 Conjunto de moléculas - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
            else:
                st.info(f"Esta funcionalidad mostrará múltiples copias de {molecula_seleccionada} en un arreglo 3D usando datos reales de ORCA.")
                
                # Verificar si existen datos de ORCA
                outputs = get_molecule_outputs(molecula_seleccionada)
                if outputs['opt'].exists():
                    st.success(f"✅ Archivo de optimización encontrado para {molecula_seleccionada}")
                    # Usar la nueva función dinámica con la molécula seleccionada
                    fig = dibujar_conjunto_molecula(molecula_seleccionada)
                    if fig:
                        st.plotly_chart(fig, width='stretch')
                else:
                    st.warning("⚠️ No se encontraron archivos de ORCA. Ejecuta primero el procesamiento.")
            
        elif option == "📦 Contenedor de moléculas":
            st.header("📦 Simulación de Contenedor Molecular")
            
            # Controles para personalizar la simulación
            col1, col2, col3 = st.columns(3)
            with col1:
                n_molecules = st.slider("Número de moléculas", 50, 1000, 400)
            with col2:
                atoms_per_molecule = st.slider("Átomos por molécula", 3, 20, 6)
            with col3:
                box_size = st.slider("Tamaño de caja (Å)", 20, 150, 80)
                
            fig = contenedor_molecular(
                n_molecules=n_molecules,
                atoms_per_molecule=atoms_per_molecule,
                box=box_size
            )
            st.plotly_chart(fig, width='stretch')

        elif option == "📈 Espectro IR":
            st.header(f"📈 Espectro Infrarrojo - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
                return
            
            fig = dibujar_espectro_ir(molecula_seleccionada)
            if fig:
                st.pyplot(fig)
        
        elif option == "⚡ Energías SCF":
            st.header(f"⚡ Análisis de Energías SCF - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
                return
            
            fig = dibujar_energias_scf(molecula_seleccionada)
            if fig:
                st.pyplot(fig)
            
        elif option == "🔬 Energías Orbitales":
            st.header(f"🔬 Análisis de Energías Orbitales - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
                return
            
            fig = dibujar_energias_orbitales(molecula_seleccionada)
            if fig:
                st.pyplot(fig)
        
        elif option == "🧬 Análisis de Población (Mulliken/Löwdin)":
            st.header(f"🧬 Análisis de Población Atómica - {molecula_seleccionada}")
            
            # Verificar que la molécula esté seleccionada
            if not molecula_seleccionada:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")
                return
            
            # Información sobre el análisis
            st.info("""
            📋 **Análisis de Población de Mulliken y Löwdin**
            
            Este análisis compara dos métodos diferentes para calcular las cargas atómicas:
            - **Mulliken**: Método tradicional que divide la densidad electrónica equitativamente entre átomos
            - **Löwdin**: Método más estable que usa orbitales ortogonalizados
            
            🔍 **¿Qué puedes analizar?**
            - Cargas atómicas de cada método
            - Diferencias entre ambos métodos
            - Distribución de población orbital
            - Estadísticas comparativas
            """)
            
            fig = dibujar_analisis_poblacion(molecula_seleccionada)
            if fig:
                st.pyplot(fig)
            
        elif option == "🔬 Trabajo de adhesión":
            st.header("🔬 Trabajo de Adhesión Molecular")
            fig = graficar_trabajo_adhesion()
            st.pyplot(fig)

        elif option == "⚛️ Molécula teórica (RDF)":
            st.header("⚛️ Función de Distribución Radial (RDF)")
            mostrar_rdf()

        elif option == "🌈 Espectros Raman":
            st.header("🌈 Análisis de Espectros IR y Raman")
            ruta_paso_3 = "modelos/FINAL_combined_spectra.txt"
            if os.path.exists(ruta_paso_3):
                mostrar_ir_raman(ruta_paso_3)
            else:
                st.error(f"No se encontró el archivo: {ruta_paso_3}")

        elif option == "🔍 Comparación con NH₃":
            st.header("🔍 Comparación Molecular vs NH₃")
            ruta = "modelos/paso_2.txt"
            if os.path.exists(ruta):
                comparar_moleculas_orca_vs_nh3(ruta)
            else:
                st.error(f"No se encontró el archivo: {ruta}")

        elif option == "🧬 IR/Raman vs NH₃":
            st.header("🧬 Comparación IR/Raman vs NH₃")
            ruta = "modelos/FINAL_combined_spectra.txt"
            if os.path.exists(ruta):
                comparar_ir_raman_vs_nh3(ruta)
            else:
                st.error(f"No se encontró el archivo: {ruta}")
        
        elif option == "📉 Desplazamientos químicos":
            st.header("📉 Análisis de Desplazamientos Químicos (RMN)")
            ruta = "modelos/paso_4.txt"
            if os.path.exists(ruta):
                comparar_rmn_s4_vs_nh3(ruta)
            else:
                st.error(f"No se encontró el archivo: {ruta}")

        elif option == "📊 Función de Distribución Radial":
            st.header("📊 Análisis de Función de Distribución Radial")
            if molecula_seleccionada:
                analyze_molecule_rdf(molecula_seleccionada)
            else:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")  

        elif option == "📉 Espectro IR Teórico":
            st.header("📉 Análisis de Espectro IR")
            if molecula_seleccionada:
                try:
                    outputs = get_molecule_outputs(molecula_seleccionada)
                    if outputs['ir-raman'].exists():
                        ir_data = parse_vibrational_frequencies(outputs['ir-raman'])
                        if ir_data is not None and not ir_data.empty:
                            # Ofrecer opciones de visualización
                            st.subheader("🔬 Opciones de Visualización")
                            
                            comparison_type = st.radio(
                                "Selecciona el tipo de análisis:",
                                ["� Comparación Educativa (Frecuencias Fundamentales vs Espectro IR)", "🧪 Solo espectro ORCA simulado"],
                                captions=[
                                    "Comparación conceptual: Modos vibracionales puros vs Simulación espectroscópica",
                                    "Muestra únicamente el espectro simulado de ORCA"
                                ]
                            )
                            
                            if comparison_type == "� Comparación Educativa (Frecuencias Fundamentales vs Espectro IR)":
                                fig = plot_ir_spectrum_comparison(ir_data, molecula_seleccionada)
                            else:
                                fig = plot_ir_spectrum(ir_data)
                                
                            st.plotly_chart(fig, width='stretch')
                            
                            # Mostrar tabla de datos
                            with st.expander("📋 Análisis Detallado de Conceptos"):
                                
                                # Información educativa
                                st.markdown("""
                                ### 📚 **Conceptos Fundamentales**
                                
                                | **Concepto** | **Qué representa** | **Dónde se ve** |
                                |--------------|-------------------|-----------------|
                                | **Frecuencias Fundamentales** | Valores numéricos exactos de modos vibracionales | Líneas verticales rojas (gráfico superior) |
                                | **Espectro IR Teórico** | Simulación de cómo se vería experimentalmente | Curva azul continua (gráfico inferior) |
                                | **Intensidad IR** | Qué tan fuerte es la absorción de cada modo | Altura de los picos en el espectro |
                                | **Modos Activos** | Vibraciones que absorben radiación infrarroja | Picos visibles en el espectro teórico |
                                """)
                                
                                col1, col2 = st.columns(2)
                                
                                with col1:
                                    st.subheader("� Frecuencias Fundamentales")
                                    st.caption("Modos vibracionales calculados por ORCA")
                                    # Crear tabla educativa
                                    freq_df = ir_data.copy()
                                    freq_df['frequency'] = freq_df['frequency'].round(1)
                                    freq_df.columns = ['Modo', 'Frecuencia (cm⁻¹)', 'Peso Relativo']
                                    freq_df['Tipo'] = 'Modo Vibracional'
                                    freq_df['Descripción'] = freq_df.apply(
                                        lambda x: f"Vibración a {x['Frecuencia (cm⁻¹)']} cm⁻¹", axis=1
                                    )
                                    st.dataframe(freq_df[['Modo', 'Frecuencia (cm⁻¹)', 'Tipo', 'Descripción']], hide_index=True)
                                
                                with col2:
                                    st.subheader("📊 Espectro IR Teórico")
                                    st.caption("Intensidades de absorción infrarroja")
                                    
                                    # Intentar leer datos de espectro IR
                                    final_ir_path = Path("modelos/FINAL_ir_spectrum.txt")
                                    if final_ir_path.exists():
                                        try:
                                            with open(final_ir_path, 'r') as f:
                                                content = f.read()
                                            
                                            ir_spec_data = []
                                            lines = content.strip().split('\n')
                                            for line in lines:
                                                if ':' in line and 'cm**-1' not in line:
                                                    parts = line.split(':')
                                                    if len(parts) >= 2:
                                                        try:
                                                            mode = int(parts[0].strip())
                                                            data_part = parts[1].strip()
                                                            values = data_part.split()
                                                            if len(values) >= 3:
                                                                freq = float(values[0])
                                                                intensity = float(values[2])
                                                                
                                                                # Clasificar intensidad
                                                                if intensity > 50:
                                                                    intensity_class = "Fuerte"
                                                                elif intensity > 10:
                                                                    intensity_class = "Media"
                                                                else:
                                                                    intensity_class = "Débil"
                                                                
                                                                ir_spec_data.append({
                                                                    'Modo': mode,
                                                                    'Frecuencia (cm⁻¹)': round(freq, 1),
                                                                    'Intensidad (km/mol)': round(intensity, 1),
                                                                    'Clasificación': intensity_class,
                                                                    'Actividad IR': 'Activo' if intensity > 1 else 'Muy débil'
                                                                })
                                                        except (ValueError, IndexError):
                                                            continue
                                            
                                            if ir_spec_data:
                                                ir_df = pd.DataFrame(ir_spec_data)
                                                st.dataframe(ir_df, hide_index=True)
                                                
                                                # Análisis comparativo educativo
                                                st.subheader("📈 Análisis Conceptual")
                                                st.markdown("""
                                                **🔍 Correspondencia entre conceptos:**
                                                - **Frecuencia:** Los valores deben coincidir entre ambas tablas
                                                - **Intensidad:** Solo visible en el espectro IR (no en frecuencias fundamentales)
                                                - **Actividad IR:** Modos con intensidad >1 km/mol son observables experimentalmente
                                                """)
                                                
                                                if len(ir_df) == len(freq_df):
                                                    # Mostrar correspondencia
                                                    corresp_df = pd.DataFrame({
                                                        'Modo': ir_df['Modo'],
                                                        'Freq. Fundamental': freq_df['Frecuencia (cm⁻¹)'].values,
                                                        'Freq. Espectro IR': ir_df['Frecuencia (cm⁻¹)'],
                                                        'Diferencia': (freq_df['Frecuencia (cm⁻¹)'].values - ir_df['Frecuencia (cm⁻¹)']).round(2),
                                                        'Intensidad IR': ir_df['Intensidad (km/mol)'],
                                                        'Observabilidad': ir_df['Actividad IR']
                                                    })
                                                    st.dataframe(corresp_df, hide_index=True)
                                                    
                                                    # Estadísticas
                                                    st.metric("Concordancia frecuencias", f"{(corresp_df['Diferencia'].abs() < 1).sum()}/{len(corresp_df)} modos")
                                                    
                                            else:
                                                st.info("No se encontraron datos de espectro IR")
                                                
                                        except Exception as e:
                                            st.error(f"Error leyendo espectro IR: {e}")
                                    else:
                                        st.info("Archivo FINAL_ir_spectrum.txt no encontrado")
                                
                        else:
                            st.warning("No se encontraron datos de frecuencias vibracionales válidos.")
                    else:
                        st.error(f"No se encontró el archivo IR/Raman para {molecula_seleccionada}")
                        st.info("💡 Ejecuta primero 'Procesar con ORCA' para generar los archivos necesarios.")
                except Exception as e:
                    st.error(f"Error procesando espectro IR: {e}")
            else:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")      
                
        elif option == "⚛️ Molécula teórica (RDF)":
            st.header("⚛️ Molécula Teórica (RDF)")
            if molecula_seleccionada:
                mostrar_rdf(molecula_seleccionada)
            else:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")

        elif option == "🌈 Espectros Raman":
            st.header("🌈 Análisis de Espectros Raman")
            if molecula_seleccionada:
                mostrar_ir_raman(molecula_seleccionada)
            else:
                st.warning("⚠️ Selecciona primero una molécula en el menú lateral.")

        elif option == "🔍 Comparación con NH₃":
            st.header("🔍 Comparación con NH₃")
            ruta = "modelos/nh3_frame.txt"
            if os.path.exists(ruta):
                comparar_moleculas_orca_vs_nh3(ruta)
            else:
                st.error(f"No se encontró el archivo: {ruta}")

        elif option == "🧬 IR/Raman vs NH₃":
            st.header("🧬 Comparación IR/Raman vs NH₃")
            ruta = "modelos/FINAL_combined_spectra.txt"
            if os.path.exists(ruta):
                comparar_ir_raman_vs_nh3(ruta)
            else:
                st.error(f"No se encontró el archivo: {ruta}")
        
        elif option == "📉 Desplazamientos químicos":
            st.header("📉 Análisis de Desplazamientos Químicos (RMN)")
            ruta = "modelos/paso_4.txt"
            if os.path.exists(ruta):
                comparar_rmn_s4_vs_nh3(ruta)
            else:
                st.error(f"No se encontró el archivo: {ruta}")

    else:
        # Si no hay molécula seleccionada
        st.warning("⚠️ Por favor, asegúrate de que la carpeta 'moleculas_xyz' contenga archivos .xyz")
        
        # Mostrar información de ayuda
        with st.expander("ℹ️ Cómo usar esta aplicación", expanded=True):
            st.write("""
            **Pasos para usar el analizador:**
            
            1. 📁 **Estructura de archivos**: Asegúrate de tener la carpeta `moleculas_xyz` con archivos .xyz
            2. 🔧 **ORCA**: Instala ORCA para ejecutar cálculos cuánticos (opcional)
            3. 🧬 **Selecciona**: Elige una molécula del menú lateral
            4. 📊 **Analiza**: Escoge el tipo de análisis que deseas realizar
            
            **Formatos soportados:**
            - Archivos XYZ estándar
            - Compatible con resultados de ORCA
            """)
            
        # Mostrar moléculas de ejemplo disponibles
        moleculas_dir = Path("moleculas_xyz")
        if moleculas_dir.exists():
            archivos_xyz = list(moleculas_dir.glob("*.xyz"))
            if archivos_xyz:
                st.subheader("🗃️ Moléculas disponibles:")
                for i, archivo in enumerate(archivos_xyz[:10]):  # Mostrar máximo 10
                    st.write(f"{i+1}. {archivo.stem}")
                if len(archivos_xyz) > 10:
                    st.write(f"... y {len(archivos_xyz) - 10} más")


if __name__ == "__main__":
    main()