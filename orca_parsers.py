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

def parse_vibrational_frequencies(out_file: Path) -> pd.DataFrame:
    """
    Extrae frecuencias vibracionales e intensidades IR de un archivo .out de ORCA.
    Si no encuentra datos en .out, intenta leer del archivo .hess correspondiente.
    Retorna DataFrame con columnas: mode, frequency, intensity
    """
    if not out_file.exists():
        return pd.DataFrame(columns=['mode', 'frequency', 'intensity'])
    
    # Primero intentar del archivo .out
    try:
        with open(out_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Extraer frecuencias del bloque VIBRATIONAL FREQUENCIES
        frequencies = {}
        freq_pattern = r"VIBRATIONAL FREQUENCIES\s*\n\s*-+\s*\n(.*?)(?=\n\s*\n|\n-+|\nNORMAL MODES|\Z)"
        freq_match = re.search(freq_pattern, content, re.DOTALL)
        
        if freq_match:
            freq_lines = freq_match.group(1).strip().split('\n')
            for line in freq_lines:
                # Formato típico: "  0:         0.00 cm**-1"
                freq_line_match = re.match(r'\s*(\d+):\s*([-+]?\d*\.?\d+)\s*cm\*\*-1', line)
                if freq_line_match:
                    mode = int(freq_line_match.group(1))
                    freq = float(freq_line_match.group(2))
                    if freq > 0:  # Solo frecuencias positivas
                        frequencies[mode] = freq
        
        # Extraer intensidades del bloque IR SPECTRUM
        intensities = {}
        ir_pattern = r"IR SPECTRUM\s*\n\s*-+.*?\n(.*?)(?=\n\s*\n|\n-+|\nRAMAN SPECTRUM|\Z)"
        ir_match = re.search(ir_pattern, content, re.DOTALL)
        
        if ir_match:
            ir_lines = ir_match.group(1).strip().split('\n')
            for line in ir_lines:
                # Formato típico: "  6:   3756.84   0.000000     42.89    0.000175"
                ir_line_match = re.match(r'\s*(\d+):\s*([-+]?\d*\.?\d+)\s+\d*\.?\d*\s+([-+]?\d*\.?\d+)', line)
                if ir_line_match:
                    mode = int(ir_line_match.group(1))
                    intensity = float(ir_line_match.group(3))
                    intensities[mode] = intensity
        
        # Combinar frecuencias e intensidades
        data = []
        for mode in frequencies:
            if mode in intensities:
                data.append({
                    'mode': mode,
                    'frequency': frequencies[mode],
                    'intensity': intensities[mode]
                })
        
        if data:
            df = pd.DataFrame(data)
            return df.sort_values('frequency').reset_index(drop=True)
        
    except Exception as e:
        print(f"Error parseando archivo .out: {e}")
    
    # Si no se encontraron datos en .out, intentar archivo .hess
    hess_file = out_file.with_suffix('.hess')
    
    if hess_file.exists():
        try:
            with open(hess_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # Buscar sección $vibrational_frequencies
            freq_pattern = r'\$vibrational_frequencies\s*\n\s*\d+\s*\n(.*?)(?=\$|$)'
            freq_match = re.search(freq_pattern, content, re.DOTALL)
            
            if freq_match:
                freq_lines = freq_match.group(1).strip().split('\n')
                
                data = []
                for line in freq_lines:
                    line = line.strip()
                    if line:
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                mode = int(parts[0])
                                freq = float(parts[1])
                                # Solo incluir frecuencias no cero (modos vibracionales reales)
                                if freq > 0.1:  # Filtrar frecuencias muy pequeñas (traslaciones/rotaciones)
                                    data.append({
                                        'mode': mode,
                                        'frequency': freq,
                                        'intensity': 1.0  # Intensidad por defecto
                                    })
                            except (ValueError, IndexError):
                                continue
                
                if data:
                    df = pd.DataFrame(data)
                    return df.sort_values('frequency').reset_index(drop=True)
                        
        except Exception as e:
            print(f"Error parseando archivo .hess: {e}")
    
    return pd.DataFrame(columns=['mode', 'frequency', 'intensity'])

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

def parse_orca_population_analysis(out_file_path: Union[str, Path]) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Extrae el análisis de población de Mulliken y Löwdin de un archivo .out de ORCA.
    Funciona para cualquier molécula.
    
    Retorna:
    {
        'mulliken': {
            'atomic_charges': {'atom_0': charge, 'atom_1': charge, ...},
            'orbital_charges': {'atom_0': {'orbital_type': charge, ...}, ...}
        },
        'loewdin': {
            'atomic_charges': {'atom_0': charge, 'atom_1': charge, ...},
            'orbital_charges': {'atom_0': {'orbital_type': charge, ...}, ...}
        }
    }
    """
    try:
        # Primero intentar leer desde archivo FINAL_population_analysis.txt
        final_file = Path("modelos/FINAL_population_analysis.txt")
        if final_file.exists():
            return parse_final_population_file(final_file)
        
        # Si no existe el archivo final, intentar extraer del .out de ORCA
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        result = {
            'mulliken': {'atomic_charges': {}, 'orbital_charges': {}},
            'loewdin': {'atomic_charges': {}, 'orbital_charges': {}}
        }
        
        # Parsear Mulliken
        result['mulliken'] = _parse_population_section(content, "MULLIKEN")
        
        # Parsear Löwdin
        result['loewdin'] = _parse_population_section(content, "LOEWDIN")
        
        return result
        
    except Exception as e:
        print(f"Error leyendo análisis de población de ORCA: {e}")
        return {
            'mulliken': {'atomic_charges': {}, 'orbital_charges': {}},
            'loewdin': {'atomic_charges': {}, 'orbital_charges': {}}
        }

def parse_final_population_file(file_path: Path) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Parsea el archivo FINAL_population_analysis.txt que contiene datos de ambos métodos.
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    result = {
        'mulliken': {'atomic_charges': {}, 'orbital_charges': {}},
        'loewdin': {'atomic_charges': {}, 'orbital_charges': {}}
    }
    
    # Dividir contenido en secciones Mulliken y Löwdin
    sections = content.split("LOEWDIN POPULATION ANALYSIS")
    
    if len(sections) >= 2:
        mulliken_content = sections[0]
        loewdin_content = "LOEWDIN POPULATION ANALYSIS" + sections[1]
        
        result['mulliken'] = _parse_population_text(mulliken_content, "MULLIKEN")
        result['loewdin'] = _parse_population_text(loewdin_content, "LOEWDIN")
    
    return result

def _parse_population_section(content: str, method: str) -> Dict[str, Dict[str, float]]:
    """
    Parsea una sección específica (MULLIKEN o LOEWDIN) del archivo .out de ORCA.
    """
    result = {'atomic_charges': {}, 'orbital_charges': {}}
    
    # Buscar sección de cargas atómicas
    atomic_pattern = rf"{method} ATOMIC CHARGES\s+-+\s+(.*?)(?=Sum of atomic charges|\n\n|\Z)"
    atomic_match = re.search(atomic_pattern, content, re.DOTALL)
    
    if atomic_match:
        atomic_lines = atomic_match.group(1).strip().split('\n')
        for line in atomic_lines:
            line = line.strip()
            if ':' in line and line.count(':') == 1:
                parts = line.split(':')
                if len(parts) == 2:
                    atom_label = parts[0].strip()
                    try:
                        charge = float(parts[1].strip())
                        result['atomic_charges'][atom_label] = charge
                    except ValueError:
                        continue
    
    # Buscar sección de cargas orbitales reducidas
    orbital_pattern = rf"{method} REDUCED ORBITAL CHARGES\s+-+\s+(.*?)(?=\n\n|\Z)"
    orbital_match = re.search(orbital_pattern, content, re.DOTALL)
    
    if orbital_match:
        orbital_lines = orbital_match.group(1).strip().split('\n')
        current_atom = None
        
        for line in orbital_lines:
            line = line.strip()
            
            # Detectar nueva sección de átomo
            if line.endswith('s:') and ' ' in line:
                current_atom = line.replace(' s:', '').strip()
                if current_atom not in result['orbital_charges']:
                    result['orbital_charges'][current_atom] = {}
                continue
            
            # Parsear orbitales
            if current_atom and ':' in line:
                parts = line.split(':')
                if len(parts) == 2:
                    orbital_type = parts[0].strip().replace(' ', '_')
                    try:
                        charge = float(parts[1].strip())
                        result['orbital_charges'][current_atom][orbital_type] = charge
                    except ValueError:
                        continue
    
    return result

def _parse_population_text(content: str, method: str) -> Dict[str, Dict[str, float]]:
    """
    Parsea el texto de población desde el archivo FINAL_population_analysis.txt.
    """
    result = {'atomic_charges': {}, 'orbital_charges': {}}
    lines = content.strip().split('\n')
    
    current_section = None
    current_atom = None
    
    for line in lines:
        line = line.strip()
        
        if 'Atomic Charges:' in line:
            current_section = 'atomic_charges'
            continue
        elif 'Reduced Orbital Charges:' in line:
            current_section = 'orbitals'
            continue
        elif line.endswith('s:') and current_section == 'orbitals':
            current_atom = line.replace(' s:', '').strip()
            if current_atom not in result['orbital_charges']:
                result['orbital_charges'][current_atom] = {}
            continue
        
        # Parsear cargas atómicas
        if current_section == 'atomic_charges' and ':' in line and not line.startswith('-'):
            parts = line.split(':')
            if len(parts) == 2:
                atom_label = parts[0].strip()
                try:
                    charge = float(parts[1].strip())
                    result['atomic_charges'][atom_label] = charge
                except ValueError:
                    continue
        
        # Parsear orbitales
        elif current_section == 'orbitals' and current_atom and ':' in line:
            parts = line.split(':')
            if len(parts) == 2:
                orbital_type = parts[0].strip().replace(' ', '_')
                try:
                    charge = float(parts[1].strip())
                    result['orbital_charges'][current_atom][orbital_type] = charge
                except ValueError:
                    continue
    
    return result

def format_population_data_for_molecule(population_data: Dict, molecule_elements: List[str]) -> Dict:
    """
    Formatea los datos de población para que sean más fáciles de usar en visualizaciones.
    Convierte etiquetas de átomos (como '0 O', '1 H') a nombres más legibles.
    """
    formatted = {
        'mulliken': {'atomic_charges': {}, 'orbital_charges': {}},
        'loewdin': {'atomic_charges': {}, 'orbital_charges': {}}
    }
    
    for method in ['mulliken', 'loewdin']:
        # Mapear cargas atómicas a nombres legibles
        element_counts = {}
        for atom_label, charge in population_data[method]['atomic_charges'].items():
            # Extraer elemento del label (ej: '0 O' -> 'O')
            element = None
            for part in atom_label.split():
                if part.isalpha():
                    element = part
                    break
            
            if element:
                if element not in element_counts:
                    element_counts[element] = 0
                
                if element_counts[element] == 0:
                    readable_name = element
                else:
                    readable_name = f"{element}{element_counts[element] + 1}"
                
                element_counts[element] += 1
                formatted[method]['atomic_charges'][readable_name] = charge
        
        # Mapear cargas orbitales
        formatted[method]['orbital_charges'] = population_data[method]['orbital_charges']
    
    return formatted