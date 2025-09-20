#!/usr/bin/env python3
"""
Script de debug para la función de parsing de energías orbitales
"""

import pandas as pd
from pathlib import Path

def parse_orca_orbital_energies(out_file_path):
    """
    Extrae las energías orbitales de un archivo .out de ORCA.
    Retorna DataFrame con datos orbitales (NO, OCC, E(Eh)).
    """
    try:
        print(f"Intentando leer desde: {out_file_path}")
        
        # Primero intentar leer desde archivo FINAL_orbital_energies.txt
        final_file = Path("modelos/FINAL_orbital_energies.txt")
        if final_file.exists():
            print(f"Encontrado archivo final: {final_file}")
            with open(final_file, 'r') as f:
                content = f.read()
            
            lines = content.strip().split('\n')
            data = []
            
            print(f"Líneas totales en archivo final: {len(lines)}")
            print("Primeras 10 líneas:")
            for i, line in enumerate(lines[:10]):
                print(f"  {i}: {repr(line)}")
            
            for line_idx, line in enumerate(lines[3:], 3):  # Saltar encabezados
                if line.strip() and not line.startswith('*'):
                    parts = line.split()
                    print(f"Línea {line_idx}: {repr(line)} -> partes: {parts}")
                    if len(parts) >= 4:
                        # Formato: NO OCC E(Eh) E(eV)
                        try:
                            orbital_data = {
                                'NO': int(parts[0]),
                                'OCC': float(parts[1]),
                                'E(Eh)': float(parts[2])
                            }
                            data.append(orbital_data)
                            print(f"  -> Agregado: {orbital_data}")
                        except (ValueError, IndexError) as e:
                            print(f"  -> Error parseando: {e}")
            
            if data:
                df = pd.DataFrame(data)
                print(f"DataFrame creado con {len(df)} orbitales")
                return df
            else:
                print("No se encontraron datos válidos en archivo final")
        
        # Si no existe el archivo final, intentar extraer del .out de ORCA
        print(f"Leyendo archivo ORCA: {out_file_path}")
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        # Buscar sección de energías orbitales
        if "ORBITAL ENERGIES" in content:
            sections = content.split("ORBITAL ENERGIES")
            print(f"Encontradas {len(sections)-1} secciones de ORBITAL ENERGIES")
            if len(sections) > 1:
                orbital_section = sections[-1]  # Tomar la última sección
                lines = orbital_section.split('\n')
                
                print(f"Líneas en sección orbital: {len(lines)}")
                print("Primeras 15 líneas de la sección:")
                for i, line in enumerate(lines[:15]):
                    print(f"  {i}: {repr(line)}")
                
                data = []
                for line_idx, line in enumerate(lines[3:], 3):  # Saltar encabezados
                    if line.strip() and not line.startswith('*'):
                        parts = line.split()
                        print(f"Línea {line_idx}: {repr(line)} -> partes: {parts}")
                        if len(parts) >= 4:
                            try:
                                # Formato típico: NO OCC E(Eh) E(eV)
                                orbital_data = {
                                    'NO': int(parts[0]),
                                    'OCC': float(parts[1]),
                                    'E(Eh)': float(parts[2])
                                }
                                data.append(orbital_data)
                                print(f"  -> Agregado: {orbital_data}")
                            except (ValueError, IndexError) as e:
                                print(f"  -> Error parseando: {e}")
                                continue
                
                if data:
                    df = pd.DataFrame(data)
                    print(f"DataFrame creado con {len(df)} orbitales desde archivo ORCA")
                    return df
                else:
                    print("No se encontraron datos válidos en archivo ORCA")
        
        # Si no encuentra datos, retornar DataFrame vacío
        print("Retornando DataFrame vacío")
        return pd.DataFrame()
        
    except Exception as e:
        print(f"Error leyendo energías orbitales de ORCA: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()

if __name__ == "__main__":
    # Probar con water
    out_path = "orca_outputs/water/water-opt.out"
    print("="*60)
    print("DEBUGGING PARSE ORBITAL ENERGIES")
    print("="*60)
    
    result = parse_orca_orbital_energies(out_path)
    print("\n" + "="*60)
    print("RESULTADO:")
    print("="*60)
    if not result.empty:
        print(f"Éxito: {len(result)} orbitales encontrados")
        print(result.head(10))
    else:
        print("Error: DataFrame vacío")
