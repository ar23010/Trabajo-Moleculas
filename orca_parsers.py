import re
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, Dict, List, Union

def parse_orca_coordinates(out_file_path: Union[str, Path]) -> Tuple[List[str], np.ndarray]:
    """
    Extrae las coordenadas finales de un archivo .out de ORCA.
    Retorna (elements, coordinates) donde coordinates es np.array de forma (n, 3)
    """
    try:
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        # Buscar la última sección de coordenadas optimizadas
        if "CARTESIAN COORDINATES (ANGSTROEM)" in content:
            # Para archivos de optimización
            sections = content.split("CARTESIAN COORDINATES (ANGSTROEM)")
            last_section = sections[-1]
        elif "CARTESIAN COORDINATES (A.U.)" in content:
            # Para otros tipos de cálculo
            sections = content.split("CARTESIAN COORDINATES (A.U.)")
            last_section = sections[-1]
        else:
            return [], np.array([])
        
        lines = last_section.split('\n')
        elements = []
        coords = []
        
        # Buscar líneas con coordenadas (formato: elemento x y z)
        for line in lines[2:]:  # Saltar header
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    element = parts[0]
                    x, y, z = map(float, parts[1:4])
                    elements.append(element)
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
            elif line.strip() == "":
                break
        
        return elements, np.array(coords) if coords else np.array([])
        
    except Exception as e:
        print(f"Error leyendo archivo ORCA: {e}")
        return [], np.array([])

def parse_orca_frequencies(out_file_path: Union[str, Path]) -> Dict[str, pd.DataFrame]:
    """
    Extrae frecuencias vibracionales e intensidades IR/Raman de un archivo .out de ORCA.
    Retorna dict con datos de IR y Raman.
    """
    try:
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        ir_data = []
        raman_data = []
        
        # Buscar sección de frecuencias
        if "VIBRATIONAL FREQUENCIES" in content:
            freq_section = content.split("VIBRATIONAL FREQUENCIES")[1].split("NORMAL MODES")[0]
            lines = freq_section.split('\n')
            
            for line in lines:
                if "cm**-1" in line and "---" not in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            freq = float(parts[1])
                            ir_data.append({'Mode': mode, 'Frequency': freq})
                        except (ValueError, IndexError):
                            continue
        
        # Buscar intensidades IR
        if "IR SPECTRUM" in content:
            ir_section = content.split("IR SPECTRUM")[1]
            if "Mode" in ir_section:
                lines = ir_section.split('\n')
                for line in lines:
                    if line.strip() and not line.startswith('Mode') and ':' in line:
                        try:
                            parts = line.replace(':', '').split()
                            if len(parts) >= 4:
                                mode = int(parts[0])
                                freq = float(parts[1])
                                intensity = float(parts[3])  # T^2 (km/mol)
                                
                                # Actualizar datos IR existentes o agregar nuevos
                                found = False
                                for item in ir_data:
                                    if item['Mode'] == mode:
                                        item['Intensity'] = intensity
                                        found = True
                                        break
                                if not found:
                                    ir_data.append({'Mode': mode, 'Frequency': freq, 'Intensity': intensity})
                        except (ValueError, IndexError):
                            continue
        
        # Buscar datos Raman
        if "RAMAN SPECTRUM" in content:
            raman_section = content.split("RAMAN SPECTRUM")[1]
            lines = raman_section.split('\n')
            for line in lines:
                if line.strip() and not line.startswith('Mode') and ':' in line:
                    try:
                        parts = line.replace(':', '').split()
                        if len(parts) >= 3:
                            mode = int(parts[0])
                            freq = float(parts[1])
                            activity = float(parts[2])
                            raman_data.append({'Mode': mode, 'Frequency': freq, 'Activity': activity})
                    except (ValueError, IndexError):
                        continue
        
        return {
            'ir': pd.DataFrame(ir_data),
            'raman': pd.DataFrame(raman_data)
        }
        
    except Exception as e:
        print(f"Error leyendo frecuencias de ORCA: {e}")
        return {'ir': pd.DataFrame(), 'raman': pd.DataFrame()}

def parse_orca_scf_energies(out_file_path: Union[str, Path]) -> Dict[str, float]:
    """
    Extrae las energías SCF de un archivo .out de ORCA.
    Retorna diccionario con las energías en Hartree y eV.
    """
    try:
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        energies = {}
        
        # Buscar sección de energías SCF
        if "TOTAL SCF ENERGY" in content:
            energy_section = content.split("TOTAL SCF ENERGY")[1].split("ORBITAL ENERGIES")[0] if "ORBITAL ENERGIES" in content else content.split("TOTAL SCF ENERGY")[1]
            
            lines = energy_section.split('\n')
            for line in lines:
                line = line.strip()
                
                # Total Energy
                if line.startswith("Total Energy"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['Total Energy'] = eh_value
                            energies['Total Energy (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # Nuclear Repulsion
                elif line.startswith("Nuclear Repulsion"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['Nuclear Repulsion'] = eh_value
                            energies['Nuclear Repulsion (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # Electronic Energy
                elif line.startswith("Electronic Energy"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['Electronic Energy'] = eh_value
                            energies['Electronic Energy (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # One Electron Energy
                elif line.startswith("One Electron Energy"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['One Electron Energy'] = eh_value
                            energies['One Electron Energy (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # Two Electron Energy
                elif line.startswith("Two Electron Energy"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['Two Electron Energy'] = eh_value
                            energies['Two Electron Energy (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # Potential Energy
                elif line.startswith("Potential Energy"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['Potential Energy'] = eh_value
                            energies['Potential Energy (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # Kinetic Energy
                elif line.startswith("Kinetic Energy"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            eh_value = float(parts[3])
                            ev_value = float(parts[5])
                            energies['Kinetic Energy'] = eh_value
                            energies['Kinetic Energy (eV)'] = ev_value
                        except (ValueError, IndexError):
                            pass
                
                # Virial Ratio
                elif line.startswith("Virial Ratio"):
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            virial_value = float(parts[2])
                            energies['Virial Ratio'] = virial_value
                        except (ValueError, IndexError):
                            pass
        
        return energies
        
    except Exception as e:
        print(f"Error leyendo energías SCF de ORCA: {e}")
        return {}

def parse_orca_orbital_energies(out_file_path: Union[str, Path]) -> pd.DataFrame:
    """
    Extrae las energías orbitales de un archivo .out de ORCA.
    Retorna DataFrame con datos orbitales (NO, OCC, E(Eh)).
    """
    try:
        # Primero intentar leer desde archivo FINAL_orbital_energies.txt
        final_file = Path("modelos/FINAL_orbital_energies.txt")
        if final_file.exists():
            with open(final_file, 'r') as f:
                content = f.read()
            
            lines = content.strip().split('\n')
            data = []
            
            for line in lines[3:]:  # Saltar encabezados
                if line.strip() and not line.startswith('*'):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            # Formato: NO OCC E(Eh) E(eV)
                            data.append({
                                'NO': int(parts[0]),
                                'OCC': float(parts[1]),
                                'E(Eh)': float(parts[2])
                            })
                        except (ValueError, IndexError):
                            # Saltar líneas que no se pueden parsear (como encabezados)
                            continue
            
            if data:
                return pd.DataFrame(data)
        
        # Si no existe el archivo final, intentar extraer del .out de ORCA
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        # Buscar sección de energías orbitales
        if "ORBITAL ENERGIES" in content:
            sections = content.split("ORBITAL ENERGIES")
            if len(sections) > 1:
                orbital_section = sections[-1]  # Tomar la última sección
                lines = orbital_section.split('\n')
                
                data = []
                for line in lines[3:]:  # Saltar encabezados
                    if line.strip() and not line.startswith('*'):
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                # Formato típico: NO OCC E(Eh) E(eV)
                                data.append({
                                    'NO': int(parts[0]),
                                    'OCC': float(parts[1]),
                                    'E(Eh)': float(parts[2])
                                })
                            except (ValueError, IndexError):
                                continue
                
                if data:
                    return pd.DataFrame(data)
        
        # Si no encuentra datos, retornar DataFrame vacío
        return pd.DataFrame()
        
    except Exception as e:
        print(f"Error leyendo energías orbitales de ORCA: {e}")
        return pd.DataFrame()

def parse_orca_xyz_file(file_path: Union[str, Path]) -> Tuple[List[str], np.ndarray]:
    """
    Parsea un archivo XYZ de ORCA.
    Retorna (elements, coordinates) donde coordinates es np.array de forma (n, 3)
    """
    elements = []
    coords = []
    in_block = False
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith("*"):
                if not in_block and line.lower().startswith("* xyz"):
                    in_block = True
                    continue
                elif in_block:
                    break
                else:
                    continue
            
            if in_block:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        element = parts[0]
                        x, y, z = map(float, parts[1:4])
                        elements.append(element)
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
    
    return elements, np.array(coords) if coords else np.array([])

def get_molecule_outputs(molecule_name: str) -> Dict[str, Path]:
    """
    Retorna los paths a los archivos .out de ORCA para una molécula específica.
    """
    base_path = Path("orca_outputs") / molecule_name
    
    return {
        'opt': base_path / f"{molecule_name}-opt.out",
        'ir-raman': base_path / f"{molecule_name}-ir-raman.out", 
        'nmr': base_path / f"{molecule_name}-nmr.out"
    }

def check_orca_outputs_exist(molecule_name: str) -> Dict[str, bool]:
    """
    Verifica qué archivos .out existen para una molécula.
    """
    outputs = get_molecule_outputs(molecule_name)
    existing = {}
    
    for calc_type, path in outputs.items():
        existing[calc_type] = path.exists()
        
    return existing