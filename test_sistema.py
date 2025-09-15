#!/usr/bin/env python3
"""
Script de prueba para verificar el funcionamiento del analizador de moléculas.
"""

import os
import sys
from pathlib import Path

def test_imports():
    """Verificar que todas las importaciones funcionan correctamente."""
    print("🔍 Verificando importaciones...")
    
    try:
        from generate_inp import OrcaInputGenerator
        print("  ✅ generate_inp importado correctamente")
    except ImportError as e:
        print(f"  ❌ Error importando generate_inp: {e}")
        return False
    
    try:
        from generar_out import OrcaOutputGenerator
        print("  ✅ generar_out importado correctamente")
    except ImportError as e:
        print(f"  ❌ Error importando generar_out: {e}")
        return False
    
    # Verificar dependencias de streamlit (sin importar streamlit para evitar conflictos)
    try:
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import plotly.graph_objects as go
        print("  ✅ Todas las dependencias disponibles")
    except ImportError as e:
        print(f"  ❌ Error con dependencias: {e}")
        return False
    
    return True

def test_directory_structure():
    """Verificar estructura de directorios."""
    print("\n📁 Verificando estructura de directorios...")
    
    required_dirs = ["moleculas_xyz"]
    created_dirs = []
    
    for dir_name in required_dirs:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            created_dirs.append(dir_name)
            print(f"  📂 Directorio creado: {dir_name}")
        else:
            print(f"  ✅ Directorio existe: {dir_name}")
    
    # Verificar archivos XYZ
    xyz_files = list(Path("moleculas_xyz").glob("*.xyz"))
    print(f"  📄 Archivos XYZ encontrados: {len(xyz_files)}")
    
    for xyz_file in xyz_files[:5]:  # Mostrar máximo 5
        print(f"    - {xyz_file.name}")
    
    if len(xyz_files) > 5:
        print(f"    ... y {len(xyz_files) - 5} más")
    
    return len(xyz_files) > 0

def test_inp_generation():
    """Probar generación de archivos .inp."""
    print("\n⚙️ Probando generación de archivos .inp...")
    
    try:
        from generate_inp import OrcaInputGenerator
        
        # Buscar un archivo XYZ para probar
        xyz_files = list(Path("moleculas_xyz").glob("*.xyz"))
        
        if not xyz_files:
            print("  ❌ No hay archivos XYZ para probar")
            return False
        
        test_xyz = xyz_files[0]
        print(f"  🧪 Probando con: {test_xyz}")
        
        # Generar archivos .inp
        generator = OrcaInputGenerator(str(test_xyz))
        generator.generate_opt_input()
        generator.generate_ir_raman_input()
        generator.generate_nmr_input()
        
        # Verificar que se crearon los archivos
        molecule_name = test_xyz.stem
        inp_dir = Path("orca_inputs") / molecule_name
        
        expected_files = [
            f"{molecule_name}-opt.inp",
            f"{molecule_name}-ir-raman.inp", 
            f"{molecule_name}-nmr.inp"
        ]
        
        success = True
        for file_name in expected_files:
            file_path = inp_dir / file_name
            if file_path.exists():
                print(f"    ✅ {file_name}")
            else:
                print(f"    ❌ {file_name} no encontrado")
                success = False
        
        return success
        
    except Exception as e:
        print(f"  ❌ Error en generación de .inp: {e}")
        return False

def test_orca_detection():
    """Verificar si ORCA está disponible."""
    print("\n🔍 Verificando disponibilidad de ORCA...")
    
    try:
        from generar_out import OrcaOutputGenerator
        
        output_gen = OrcaOutputGenerator()
        
        if output_gen.orca_path:
            print(f"  ✅ ORCA encontrado en: {output_gen.orca_path}")
            return True
        else:
            print("  ⚠️  ORCA no encontrado. Los cálculos no se pueden ejecutar automáticamente.")
            print("     Los archivos .inp se pueden ejecutar manualmente.")
            return False
            
    except Exception as e:
        print(f"  ❌ Error verificando ORCA: {e}")
        return False

def test_file_parsing():
    """Probar parseo de archivos XYZ."""
    print("\n📖 Probando parseo de archivos XYZ...")
    
    try:
        xyz_files = list(Path("moleculas_xyz").glob("*.xyz"))
        
        if not xyz_files:
            print("  ❌ No hay archivos XYZ para probar")
            return False
            
        test_xyz = xyz_files[0]
        print(f"  📄 Probando con: {test_xyz}")
        
        with open(test_xyz, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            print("  ❌ Archivo XYZ inválido (muy pocas líneas)")
            return False
        
        try:
            num_atoms = int(lines[0].strip())
            description = lines[1].strip()
            print(f"    ✅ Átomos: {num_atoms}")
            print(f"    ✅ Descripción: {description}")
            
            # Parsear algunas coordenadas
            atoms_parsed = 0
            for i in range(2, min(len(lines), num_atoms + 2)):
                parts = lines[i].strip().split()
                if len(parts) >= 4:
                    element = parts[0]
                    coords = list(map(float, parts[1:4]))
                    atoms_parsed += 1
            
            print(f"    ✅ Coordenadas parseadas: {atoms_parsed}/{num_atoms}")
            return atoms_parsed > 0
            
        except ValueError as e:
            print(f"  ❌ Error parseando archivo: {e}")
            return False
            
    except Exception as e:
        print(f"  ❌ Error en parseo: {e}")
        return False

def run_all_tests():
    """Ejecutar todas las pruebas."""
    print("🧪 INICIANDO PRUEBAS DEL ANALIZADOR DE MOLÉCULAS")
    print("=" * 50)
    
    tests = [
        ("Importaciones", test_imports),
        ("Estructura de directorios", test_directory_structure),
        ("Generación de archivos .inp", test_inp_generation),
        ("Detección de ORCA", test_orca_detection),
        ("Parseo de archivos", test_file_parsing),
    ]
    
    results = {}
    
    for test_name, test_func in tests:
        try:
            results[test_name] = test_func()
        except Exception as e:
            print(f"\n❌ Error inesperado en {test_name}: {e}")
            results[test_name] = False
    
    # Resumen final
    print("\n" + "=" * 50)
    print("📊 RESUMEN DE PRUEBAS")
    print("=" * 50)
    
    passed = 0
    total = len(tests)
    
    for test_name, result in results.items():
        status = "✅ PASÓ" if result else "❌ FALLÓ"
        print(f"{status:10} {test_name}")
        if result:
            passed += 1
    
    print("-" * 50)
    print(f"Pruebas pasadas: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 ¡Todas las pruebas pasaron! El sistema está listo.")
        print("\nPuedes ejecutar la aplicación con:")
        print("streamlit run spectra_app.py")
    else:
        print(f"\n⚠️  {total - passed} pruebas fallaron. Revisa los errores arriba.")
    
    return passed == total

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
