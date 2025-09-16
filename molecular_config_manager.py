"""
Sistema de configuración automática para moléculas.
Este módulo maneja la actualización automática de los archivos de pasos
basándose en la molécula seleccionada.
"""
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np

class MolecularConfigManager:
    """Gestor de configuración automática para moléculas."""
    
    def __init__(self, base_dir: str = "."):
        """Inicializar el gestor de configuración."""
        self.base_dir = Path(base_dir)
        self.moleculas_dir = self.base_dir / "moleculas_xyz"
        self.modelos_dir = self.base_dir / "modelos"
        self.orca_inputs_dir = self.base_dir / "orca_inputs"
        
        # Archivos de pasos que necesitan ser actualizados
        self.paso_files = {
            1: self.modelos_dir / "paso_1.txt",
            2: self.modelos_dir / "paso_2.txt", 
            4: self.modelos_dir / "paso_4.txt"
        }
        
        # Verificar que existen los directorios necesarios
        self._verify_directories()
        
    def _verify_directories(self):
        """Verificar que existen los directorios necesarios."""
        for directory in [self.moleculas_dir, self.modelos_dir]:
            if not directory.exists():
                raise FileNotFoundError(f"Directorio requerido no encontrado: {directory}")
    
    def get_available_molecules(self) -> List[str]:
        """Obtener lista de moléculas disponibles."""
        if not self.moleculas_dir.exists():
            return []
        
        molecules = []
        for xyz_file in self.moleculas_dir.glob("*.xyz"):
            molecules.append(xyz_file.stem)
        
        return sorted(molecules)
    
    def read_xyz_file(self, molecule_name: str) -> Tuple[List[str], np.ndarray, str]:
        """
        Leer archivo XYZ y extraer información.
        
        Returns:
            Tuple[elementos, coordenadas, descripción]
        """
        xyz_path = self.moleculas_dir / f"{molecule_name}.xyz"
        
        if not xyz_path.exists():
            raise FileNotFoundError(f"Archivo XYZ no encontrado: {xyz_path}")
        
        elements = []
        coordinates = []
        description = ""
        
        with open(xyz_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            raise ValueError(f"Archivo XYZ inválido: {xyz_path}")
        
        # Primera línea: número de átomos
        num_atoms = int(lines[0].strip())
        
        # Segunda línea: descripción
        if len(lines) > 1:
            description = lines[1].strip()
        
        # Líneas siguientes: elementos y coordenadas
        for i in range(2, min(len(lines), num_atoms + 2)):
            parts = lines[i].strip().split()
            if len(parts) >= 4:
                element = parts[0]
                x, y, z = map(float, parts[1:4])
                elements.append(element)
                coordinates.append([x, y, z])
        
        return elements, np.array(coordinates), description
    
    def generate_xyz_block(self, elements: List[str], coordinates: np.ndarray, 
                          charge: int = 0, multiplicity: int = 1) -> str:
        """Generar bloque XYZ para archivos ORCA."""
        xyz_block = f"* xyz {charge} {multiplicity}\n"
        
        for element, coord in zip(elements, coordinates):
            xyz_block += f"{element:2s}   {coord[0]:12.6f}   {coord[1]:12.6f}   {coord[2]:12.6f}\n"
        
        xyz_block += "*\n"
        return xyz_block
    
    def update_paso_1(self, molecule_name: str, elements: List[str], 
                      coordinates: np.ndarray) -> bool:
        """
        Actualizar paso_1.txt con los datos de la molécula.
        Este archivo contiene configuración de optimización de geometría.
        """
        try:
            paso_1_path = self.paso_files[1]
            
            # Plantilla para paso 1 (optimización)
            template = """!B3LYP 6-31+G(d,p) OPT DEFGRID2 TIGHTSCF D3BJ PAL8
!geom
  iter 200
end
{xyz_block}
"""
            
            # Generar bloque XYZ
            xyz_block = self.generate_xyz_block(elements, coordinates)
            
            # Escribir archivo actualizado
            content = template.format(xyz_block=xyz_block)
            
            with open(paso_1_path, 'w') as f:
                f.write(content)
            
            print(f"✅ Actualizado paso_1.txt con datos de {molecule_name}")
            return True
            
        except Exception as e:
            print(f"❌ Error actualizando paso_1.txt: {e}")
            return False
    
    def update_paso_2(self, molecule_name: str, elements: List[str], 
                      coordinates: np.ndarray) -> bool:
        """
        Actualizar paso_2.txt con los datos de la molécula.
        Este archivo contiene configuración para espectros IR/Raman.
        """
        try:
            paso_2_path = self.paso_files[2]
            
            # Plantilla para paso 2 (IR/Raman)
            template = """# avogadro generated ORCA input file
# Basic Mode
#
! B3LYP 6-31+G(d,p) NumFreq DEFGRID2 TIGHTSCF D3BJ PAL8
%elprop Polar 1 end
{xyz_block}
"""
            
            # Generar bloque XYZ
            xyz_block = self.generate_xyz_block(elements, coordinates)
            
            # Escribir archivo actualizado
            content = template.format(xyz_block=xyz_block)
            
            with open(paso_2_path, 'w') as f:
                f.write(content)
            
            print(f"✅ Actualizado paso_2.txt con datos de {molecule_name}")
            return True
            
        except Exception as e:
            print(f"❌ Error actualizando paso_2.txt: {e}")
            return False
    
    def update_paso_4(self, molecule_name: str, elements: List[str], 
                      coordinates: np.ndarray) -> bool:
        """
        Actualizar paso_4.txt con los datos de la molécula.
        Este archivo contiene configuración para RMN.
        """
        try:
            paso_4_path = self.paso_files[4]
            
            # Plantilla para paso 4 (RMN)
            template = """! B3LYP 6-311+G(2d,p) AutoAux DEFGRID2 TIGHTSCF D3BJ
%pal nprocs 8 end
# Opcional: NMR directo en lugar de EPRNMR
# ! NMR

{xyz_block}
%EPRNMR
  # Calcula blindajes para todos los núcleos (incluye 1H y 13C).
  # Si quisieras limitar: NUCLEI = H, C
  NUCLEI = ALL {{sSall}}
END

"""
            
            # Generar bloque XYZ
            xyz_block = self.generate_xyz_block(elements, coordinates)
            
            # Escribir archivo actualizado
            content = template.format(xyz_block=xyz_block)
            
            with open(paso_4_path, 'w') as f:
                f.write(content)
            
            print(f"✅ Actualizado paso_4.txt con datos de {molecule_name}")
            return True
            
        except Exception as e:
            print(f"❌ Error actualizando paso_4.txt: {e}")
            return False
    
    def auto_configure_molecule(self, molecule_name: str) -> Dict[str, bool]:
        """
        Configurar automáticamente todos los archivos de pasos para una molécula.
        
        Args:
            molecule_name: Nombre de la molécula (sin extensión .xyz)
        
        Returns:
            Dict con el resultado de cada actualización
        """
        results = {}
        
        try:
            print(f"🔄 Configurando automáticamente: {molecule_name}")
            
            # Leer datos de la molécula
            elements, coordinates, description = self.read_xyz_file(molecule_name)
            
            print(f"📊 Molécula: {molecule_name}")
            print(f"📄 Descripción: {description}")
            print(f"⚛️ Átomos: {len(elements)} ({', '.join(set(elements))})")
            
            # Actualizar cada paso
            results['paso_1'] = self.update_paso_1(molecule_name, elements, coordinates)
            results['paso_2'] = self.update_paso_2(molecule_name, elements, coordinates)
            results['paso_4'] = self.update_paso_4(molecule_name, elements, coordinates)
            
            # Resumen
            successful_updates = sum(results.values())
            total_updates = len(results)
            
            print(f"\n📈 Resumen de configuración:")
            print(f"   ✅ Exitosos: {successful_updates}/{total_updates}")
            
            if successful_updates == total_updates:
                print(f"🎉 Configuración completa para {molecule_name}")
            else:
                print(f"⚠️ Configuración parcial para {molecule_name}")
                
        except Exception as e:
            print(f"❌ Error en configuración automática: {e}")
            results = {'paso_1': False, 'paso_2': False, 'paso_4': False}
        
        return results
    
    def create_molecule_config_file(self, molecule_name: str) -> bool:
        """Crear archivo de configuración específico para una molécula."""
        try:
            config_path = self.modelos_dir / f"config_{molecule_name}.txt"
            
            elements, coordinates, description = self.read_xyz_file(molecule_name)
            
            # Información de configuración
            config_content = f"""# Configuración automática para {molecule_name}
# Generado automáticamente por MolecularConfigManager
# Fecha: {os.popen('date').read().strip()}

MOLECULE_NAME: {molecule_name}
DESCRIPTION: {description}
NUM_ATOMS: {len(elements)}
ELEMENTS: {', '.join(set(elements))}

# Coordenadas (Angstrom)
COORDINATES:
"""
            for i, (element, coord) in enumerate(zip(elements, coordinates)):
                config_content += f"{i+1:3d}  {element:2s}   {coord[0]:12.6f}   {coord[1]:12.6f}   {coord[2]:12.6f}\n"
            
            with open(config_path, 'w') as f:
                f.write(config_content)
            
            print(f"✅ Archivo de configuración creado: {config_path}")
            return True
            
        except Exception as e:
            print(f"❌ Error creando archivo de configuración: {e}")
            return False
    
    def get_current_molecule_from_paso(self, paso_num: int = 2) -> Optional[str]:
        """
        Intentar determinar qué molécula está configurada actualmente
        analizando el archivo de paso.
        """
        try:
            paso_path = self.paso_files[paso_num]
            
            if not paso_path.exists():
                return None
                
            with open(paso_path, 'r') as f:
                content = f.read()
            
            # Extraer coordenadas del archivo
            lines = content.split('\n')
            current_elements = []
            current_coords = []
            
            in_xyz_block = False
            for line in lines:
                line = line.strip()
                
                if line.startswith('* xyz'):
                    in_xyz_block = True
                    continue
                elif line.startswith('*') and in_xyz_block:
                    break
                elif in_xyz_block and line:
                    parts = line.split()
                    if len(parts) >= 4:
                        element = parts[0]
                        coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                        current_elements.append(element)
                        current_coords.append(coords)
            
            # Comparar con moléculas disponibles
            current_coords = np.array(current_coords)
            
            for molecule_name in self.get_available_molecules():
                try:
                    mol_elements, mol_coords, _ = self.read_xyz_file(molecule_name)
                    
                    # Comparar número de átomos y tipos
                    if len(mol_elements) == len(current_elements):
                        if sorted(mol_elements) == sorted(current_elements):
                            # Comparar coordenadas (con tolerancia)
                            if mol_coords.shape == current_coords.shape:
                                diff = np.abs(mol_coords - current_coords)
                                if np.all(diff < 0.001):  # Tolerancia de 0.001 Å
                                    return molecule_name
                                    
                except Exception:
                    continue
            
            return None
            
        except Exception as e:
            print(f"❌ Error determinando molécula actual: {e}")
            return None


def main():
    """Función principal para testing del módulo."""
    print("🧬 Gestor de Configuración Molecular")
    print("=" * 50)
    
    try:
        manager = MolecularConfigManager()
        
        # Listar moléculas disponibles
        molecules = manager.get_available_molecules()
        print(f"📋 Moléculas disponibles: {len(molecules)}")
        for i, mol in enumerate(molecules, 1):
            print(f"   {i}. {mol}")
        
        # Ejemplo de uso
        if molecules:
            example_molecule = molecules[0]
            print(f"\n🔬 Ejemplo con {example_molecule}:")
            
            # Configurar automáticamente
            results = manager.auto_configure_molecule(example_molecule)
            
            # Crear archivo de configuración
            manager.create_molecule_config_file(example_molecule)
            
            # Detectar molécula actual
            current = manager.get_current_molecule_from_paso(2)
            if current:
                print(f"\n🎯 Molécula actualmente configurada: {current}")
            else:
                print(f"\n❓ No se pudo determinar la molécula actual")
        
    except Exception as e:
        print(f"❌ Error en el sistema: {e}")


if __name__ == "__main__":
    main()
