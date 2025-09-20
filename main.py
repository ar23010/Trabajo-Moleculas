# main.py
import streamlit as st
from pathlib import Path
import os

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
            "🔬 Trabajo de adhesión",
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
                st.plotly_chart(fig, use_container_width=True)
            
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
                        st.plotly_chart(fig, use_container_width=True)
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
            st.plotly_chart(fig, use_container_width=True)

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