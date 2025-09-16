#!/usr/bin/env python3
"""
Script de demostración del sistema de configuración automática.
Prueba la funcionalidad del MolecularConfigManager.
"""

import os
import sys
from pathlib import Path

# Agregar el directorio actual al path
sys.path.append(str(Path(__file__).parent))

from molecular_config_manager import MolecularConfigManager

def demo_configuracion_automatica():
    """Demostración del sistema de configuración automática."""
    print("🧬 DEMOSTRACIÓN - Sistema de Configuración Automática")
    print("=" * 60)
    
    try:
        # Inicializar gestor
        manager = MolecularConfigManager()
        print("✅ Gestor de configuración inicializado")
        
        # Listar moléculas disponibles
        molecules = manager.get_available_molecules()
        print(f"\n📋 Moléculas disponibles: {len(molecules)}")
        for i, mol in enumerate(molecules, 1):
            print(f"   {i}. {mol}")
        
        if not molecules:
            print("❌ No se encontraron moléculas en moleculas_xyz/")
            return
        
        # Detectar molécula actual
        current_molecule = manager.get_current_molecule_from_paso(2)
        if current_molecule:
            print(f"\n🎯 Molécula actualmente configurada: {current_molecule}")
        else:
            print(f"\n❓ No se detectó ninguna molécula configurada")
        
        # Seleccionar molécula para demostración
        test_molecule = molecules[0]  # Usar la primera molécula
        print(f"\n🔬 Probando con: {test_molecule}")
        
        # Mostrar información de la molécula
        elements, coordinates, description = manager.read_xyz_file(test_molecule)
        print(f"   📊 Átomos: {len(elements)}")
        print(f"   📄 Descripción: {description}")
        print(f"   ⚛️ Elementos: {', '.join(set(elements))}")
        
        # Estado inicial de archivos
        print(f"\n📁 Estado inicial de archivos:")
        for paso, path in manager.paso_files.items():
            exists = "✅" if path.exists() else "❌"
            print(f"   paso_{paso}.txt: {exists}")
        
        # Configuración automática
        print(f"\n🔄 Configurando automáticamente {test_molecule}...")
        results = manager.auto_configure_molecule(test_molecule)
        
        # Resultados
        print(f"\n📈 Resultados:")
        successful = sum(results.values())
        total = len(results)
        
        for paso, success in results.items():
            status = "✅" if success else "❌"
            print(f"   {paso}: {status}")
        
        print(f"\n🎉 Configuración completada: {successful}/{total} archivos")
        
        # Estado final
        print(f"\n📁 Estado final de archivos:")
        for paso, path in manager.paso_files.items():
            if path.exists():
                size = path.stat().st_size
                print(f"   paso_{paso}.txt: ✅ ({size} bytes)")
            else:
                print(f"   paso_{paso}.txt: ❌")
        
        # Verificar detección
        print(f"\n🔍 Verificando detección automática...")
        detected = manager.get_current_molecule_from_paso(2)
        if detected == test_molecule:
            print(f"   ✅ Detección correcta: {detected}")
        else:
            print(f"   ⚠️ Detección inconsistente: {detected} vs {test_molecule}")
        
        # Crear archivo de configuración
        print(f"\n📄 Creando archivo de configuración...")
        config_created = manager.create_molecule_config_file(test_molecule)
        if config_created:
            print(f"   ✅ Archivo de configuración creado")
        else:
            print(f"   ❌ Error creando archivo de configuración")
            
        print(f"\n🎯 Demostración completada exitosamente")
        
    except Exception as e:
        print(f"❌ Error en la demostración: {e}")
        import traceback
        traceback.print_exc()

def demo_reactividad():
    """Demostración de la reactividad del sistema."""
    print("\n" + "=" * 60)
    print("🔄 DEMOSTRACIÓN - Reactividad del Sistema")
    print("=" * 60)
    
    try:
        manager = MolecularConfigManager()
        molecules = manager.get_available_molecules()
        
        if len(molecules) < 2:
            print("❌ Se necesitan al menos 2 moléculas para probar reactividad")
            return
        
        # Configurar primera molécula
        mol1 = molecules[0]
        print(f"🧪 Configurando primera molécula: {mol1}")
        manager.auto_configure_molecule(mol1)
        
        current = manager.get_current_molecule_from_paso(2)
        print(f"✅ Molécula detectada: {current}")
        
        # Cambiar a segunda molécula
        mol2 = molecules[1]
        print(f"\n🔄 Cambiando a segunda molécula: {mol2}")
        manager.auto_configure_molecule(mol2)
        
        current = manager.get_current_molecule_from_paso(2)
        print(f"✅ Nueva molécula detectada: {current}")
        
        # Verificar que el cambio fue detectado correctamente
        if current == mol2:
            print(f"🎉 Reactividad funcionando correctamente!")
        else:
            print(f"⚠️ Problema con la reactividad")
        
    except Exception as e:
        print(f"❌ Error en demostración de reactividad: {e}")

def verificar_integracion_streamlit():
    """Verificar que la integración con Streamlit funcione."""
    print("\n" + "=" * 60)
    print("🖥️ VERIFICACIÓN - Integración con Streamlit")
    print("=" * 60)
    
    try:
        # Verificar imports
        import streamlit as st
        print("✅ Streamlit disponible")
        
        # Verificar que spectra_app puede importar el nuevo módulo
        try:
            from spectra_app import seleccionar_molecula
            print("✅ Función reactiva de selección disponible")
        except ImportError as e:
            print(f"❌ Error importando función reactiva: {e}")
        
        # Verificar archivos necesarios
        required_files = [
            "spectra_app.py",
            "molecular_config_manager.py",
            "generate_inp.py", 
            "generar_out.py"
        ]
        
        for file in required_files:
            if Path(file).exists():
                print(f"✅ {file} encontrado")
            else:
                print(f"❌ {file} no encontrado")
        
        # Verificar directorios
        required_dirs = [
            "moleculas_xyz",
            "modelos", 
            "orca_inputs",
            "orca_outputs"
        ]
        
        for dir in required_dirs:
            if Path(dir).exists():
                print(f"✅ Directorio {dir}/ encontrado")
            else:
                print(f"⚠️ Directorio {dir}/ no encontrado (se creará automáticamente)")
        
        print("🎯 Verificación de integración completada")
        
    except Exception as e:
        print(f"❌ Error en verificación: {e}")

def main():
    """Función principal del script de demostración."""
    print("🚀 Iniciando demostración del sistema de configuración automática...")
    
    # Cambiar al directorio del script
    script_dir = Path(__file__).parent
    os.chdir(script_dir)
    
    # Ejecutar demostraciones
    demo_configuracion_automatica()
    demo_reactividad()
    verificar_integracion_streamlit()
    
    print("\n" + "=" * 60)
    print("📝 INSTRUCCIONES PARA USAR EL SISTEMA:")
    print("=" * 60)
    print("1. Ejecuta la aplicación Streamlit:")
    print("   streamlit run spectra_app.py")
    print("")
    print("2. En el sidebar, selecciona una molécula")
    print("   → El sistema automáticamente actualizará los archivos de pasos")
    print("")
    print("3. Ve al 'Estado de configuración' para ver el estado actual")
    print("")
    print("4. Los archivos paso_1.txt, paso_2.txt y paso_4.txt se actualizan")
    print("   automáticamente cada vez que cambias de molécula")
    print("")
    print("🎉 ¡El sistema está listo para usar!")

if __name__ == "__main__":
    main()
