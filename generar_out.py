import os
import sys
import subprocess
import shutil
import re
from pathlib import Path
from typing import List, Optional, Dict, Any
import time

# Ruta del ejecutable de ORCA (personalizable por el usuario)
ORCA_BIN = "/home/fercho/orca_6_1_0_linux_x86-64_shared_openmpi418/orca"

class OrcaSpectrumExtractor:
    """
    Extrae información de espectros IR y Raman de archivos .out de ORCA
    """
    
    @staticmethod
    def extract_ir_spectrum(out_file_content: str) -> List[Dict[str, Any]]:
        """Extrae la tabla de espectro IR del contenido del archivo .out"""
        ir_spectrum = []
        
        # Patrón para encontrar la sección IR SPECTRUM
        ir_pattern = r"IR SPECTRUM\s+-+-\s+(Mode\s+freq\s+eps\s+Int\s+T\*\*2.*?)(?:\n\n|\Z)"
        ir_match = re.search(ir_pattern, out_file_content, re.DOTALL)
        
        if not ir_match:
            return ir_spectrum
            
        ir_table = ir_match.group(1)
        
        # Extraer líneas de datos (saltando encabezados)
        lines = ir_table.split('\n')
        data_started = False
        
        for line in lines:
            if line.startswith('---'):
                data_started = True
                continue
                
            if data_started and line.strip() and not line.startswith('Mode'):
                # Parsear línea de datos IR
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        mode = parts[0].rstrip(':')
                        freq = float(parts[1])
                        eps = float(parts[2])
                        intensity = float(parts[3])
                        t_squared = float(parts[4])
                        
                        # Extraer componentes del vector (pueden estar entre paréntesis)
                        tx, ty, tz = 0.0, 0.0, 0.0
                        if '(' in line and ')' in line:
                            vector_part = re.search(r'\((.*?)\)', line)
                            if vector_part:
                                vector_parts = vector_part.group(1).split()
                                if len(vector_parts) >= 3:
                                    tx = float(vector_parts[0])
                                    ty = float(vector_parts[1])
                                    tz = float(vector_parts[2])
                        
                        ir_spectrum.append({
                            'mode': mode,
                            'frequency': freq,
                            'epsilon': eps,
                            'intensity': intensity,
                            't_squared': t_squared,
                            'tx': tx,
                            'ty': ty,
                            'tz': tz
                        })
                    except (ValueError, IndexError):
                        continue
        
        return ir_spectrum
    
    @staticmethod
    def extract_raman_spectrum(out_file_content: str) -> List[Dict[str, Any]]:
        """Extrae la tabla de espectro Raman del contenido del archivo .out"""
        raman_spectrum = []
        
        # Patrón para encontrar la sección RAMAN SPECTRUM
        raman_pattern = r"RAMAN SPECTRUM\s+-+-\s+(Mode\s+freq.*?Activity.*?Depolarization.*?)(?:\n\n|\Z)"
        raman_match = re.search(raman_pattern, out_file_content, re.DOTALL)
        
        if not raman_match:
            return raman_spectrum
            
        raman_table = raman_match.group(1)
        
        # Extraer líneas de datos
        lines = raman_table.split('\n')
        data_started = False
        
        for line in lines:
            if line.startswith('---'):
                data_started = True
                continue
                
            if data_started and line.strip() and not line.startswith('Mode'):
                # Parsear línea de datos Raman
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        mode = parts[0].rstrip(':')
                        freq = float(parts[2])
                        activity = float(parts[3])
                        depolarization = float(parts[4])
                        
                        raman_spectrum.append({
                            'mode': mode,
                            'frequency': freq,
                            'activity': activity,
                            'depolarization': depolarization
                        })
                    except (ValueError, IndexError):
                        continue
        
        return raman_spectrum
    
    @staticmethod
    def extract_spectra_from_file(out_file_path: Path) -> Dict[str, Any]:
        """Extrae espectros IR y Raman de un archivo .out específico"""
        try:
            with open(out_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            return {
                'ir_spectrum': OrcaSpectrumExtractor.extract_ir_spectrum(content),
                'raman_spectrum': OrcaSpectrumExtractor.extract_raman_spectrum(content)
            }
        except Exception as e:
            print(f"❌ Error leyendo archivo {out_file_path}: {e}")
            return {'ir_spectrum': [], 'raman_spectrum': []}


class OrcaOutputGenerator:
    """
    Generador de archivos .out de ORCA a partir de archivos .inp existentes.
    Ejecuta cálculos ORCA secuencialmente: Optimización, IR/Raman, RMN.
    """

    def __init__(self, base_inp_dir: str = "orca_inputs", base_out_dir: str = "orca_outputs"):
        self.base_inp_dir = Path(base_inp_dir)
        self.base_out_dir = Path(base_out_dir)
        self.spectrum_extractor = OrcaSpectrumExtractor()
        
        # Verificar que ORCA esté disponible
        self.orca_path = self._find_orca()
        
        # Crear carpeta de salida
        self.base_out_dir.mkdir(exist_ok=True)
        

    def _find_orca(self) -> Optional[str]:
        """Busca el ejecutable de ORCA en el sistema o usa ruta personalizada."""
        # Si el usuario define ORCA_BIN, usar esa ruta
        if ORCA_BIN:
            if os.path.isfile(ORCA_BIN) and os.access(ORCA_BIN, os.X_OK):
                print(f"✅ ORCA personalizado encontrado en: {ORCA_BIN}")
                return ORCA_BIN
            else:
                print(f"❌ Ruta personalizada de ORCA no válida: {ORCA_BIN}")
                return None
        # Si no, buscar en ubicaciones comunes y PATH
        possible_paths = [
            "/opt/orca/orca",
            "/usr/local/bin/orca",
            "/usr/bin/orca",
            "orca"  # Si está en PATH
        ]
        for path in possible_paths:
            found = shutil.which(path)
            if found:
                print(f"✅ ORCA encontrado en: {found}")
                return found
        print("⚠️  ORCA no encontrado. Asegúrate de que esté instalado y en PATH o define ORCA_BIN.")
        return None
        
    def _run_orca_calculation(self, inp_file: Path, out_dir: Path) -> bool:
        """Ejecuta un cálculo ORCA individual."""
        if not self.orca_path:
            print("❌ ORCA no está disponible")
            return False
            
        try:
            # Crear carpeta de salida específica
            out_dir.mkdir(exist_ok=True)
            
            # Nombres de archivos
            inp_name = inp_file.stem
            out_file = out_dir / f"{inp_name}.out"
            
            # Copiar archivo .inp a la carpeta de salida
            inp_copy = out_dir / inp_file.name
            shutil.copy2(inp_file, inp_copy)
            
            print(f"🔄 Ejecutando cálculo ORCA: {inp_file.name}")
            
            # Ejecutar ORCA con la ruta absoluta del archivo
            cmd = [self.orca_path, inp_copy.name]  # Solo el nombre del archivo
            
            # Mostrar el comando que se va a ejecutar (debug)
            print(f"📍 Ejecutando comando: {' '.join(cmd)} en directorio: {out_dir}")
            
            with open(out_file, "w") as f:
                process = subprocess.run(
                    cmd,
                    cwd=str(out_dir),  # Convertir Path a string y ejecutar desde el directorio de salida
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600  # 1 hora máximo por cálculo
                )
            
            if process.returncode == 0:
                print(f"✅ Cálculo completado: {out_file}")
                
                # Extraer espectros después del cálculo exitoso
                spectra = self.spectrum_extractor.extract_spectra_from_file(out_file)
                self._save_spectra_data(spectra, out_dir, inp_name)
                
                return True
            else:
                print(f"❌ Error en cálculo: {inp_file.name}")
                print(f"Error: {process.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            print(f"⏰ Timeout en cálculo: {inp_file.name}")
            return False
        except Exception as e:
            print(f"❌ Error inesperado: {e}")
            return False
    
    def _save_spectra_data(self, spectra: Dict[str, Any], out_dir: Path, base_name: str):
        """Guarda los datos de espectros en archivos separados"""
        # Guardar datos IR
        if spectra['ir_spectrum']:
            ir_file = out_dir / f"{base_name}_ir_spectrum.txt"
            with open(ir_file, 'w') as f:
                f.write("IR SPECTRUM\n")
                f.write("Mode\tFrequency (cm⁻¹)\tEpsilon\tIntensity (km/mol)\tT²\n")
                for peak in spectra['ir_spectrum']:
                    f.write(f"{peak['mode']}\t{peak['frequency']:.2f}\t"
                           f"{peak['epsilon']:.6f}\t{peak['intensity']:.2f}\t"
                           f"{peak['t_squared']:.6f}\n")
            print(f"📊 Datos IR guardados en: {ir_file}")
        
        # Guardar datos Raman
        if spectra['raman_spectrum']:
            raman_file = out_dir / f"{base_name}_raman_spectrum.txt"
            with open(raman_file, 'w') as f:
                f.write("RAMAN SPECTRUM\n")
                f.write("Mode\tFrequency (cm⁻¹)\tActivity\tDepolarization\n")
                for peak in spectra['raman_spectrum']:
                    f.write(f"{peak['mode']}\t{peak['frequency']:.2f}\t"
                           f"{peak['activity']:.6f}\t{peak['depolarization']:.6f}\n")
            print(f"📊 Datos Raman guardados en: {raman_file}")
            
    def process_molecule(self, molecule_name: str) -> dict:
        """
        Procesa todos los cálculos ORCA para una molécula específica.
        
        Args:
            molecule_name: Nombre de la molécula (sin extensión)
            
        Returns:
            dict: Resultados de los cálculos
        """
        molecule_inp_dir = self.base_inp_dir / molecule_name
        molecule_out_dir = self.base_out_dir / molecule_name
        
        if not molecule_inp_dir.exists():
            print(f"❌ No se encuentra carpeta de inputs: {molecule_inp_dir}")
            return {"success": False, "error": "Inputs no encontrados"}
            
        # Crear carpeta de salida para la molécula
        molecule_out_dir.mkdir(exist_ok=True)
        
        results = {
            "molecule": molecule_name,
            "success": True,
            "calculations": {},
            "spectra": {}
        }
        
        # Orden de cálculos: opt -> ir-raman -> nmr
        calc_order = ["opt", "ir-raman", "nmr"]
        
        for calc_type in calc_order:
            inp_pattern = f"{molecule_name}-{calc_type}.inp"
            inp_file = molecule_inp_dir / inp_pattern
            
            if inp_file.exists():
                print(f"\n📊 Iniciando cálculo {calc_type.upper()} para {molecule_name}")
                start_time = time.time()
                
                success = self._run_orca_calculation(inp_file, molecule_out_dir)
                
                elapsed_time = time.time() - start_time
                results["calculations"][calc_type] = {
                    "success": success,
                    "time": elapsed_time,
                    "output_file": str(molecule_out_dir / f"{molecule_name}-{calc_type}.out")
                }
                
                # Extraer espectros solo para cálculos IR-Raman exitosos
                if calc_type == "ir-raman" and success:
                    out_file = molecule_out_dir / f"{molecule_name}-{calc_type}.out"
                    spectra = self.spectrum_extractor.extract_spectra_from_file(out_file)
                    results["spectra"] = spectra
                
                if not success:
                    results["success"] = False
                    print(f"⚠️  Continuando con siguiente cálculo...")
                    
            else:
                print(f"⚠️  Archivo no encontrado: {inp_file}")
                results["calculations"][calc_type] = {
                    "success": False,
                    "error": "Archivo input no encontrado"
                }
                
        return results
        
    def process_all_molecules(self) -> dict:
        """Procesa todos los cálculos ORCA para todas las moléculas disponibles."""
        if not self.base_inp_dir.exists():
            return {"success": False, "error": "Carpeta orca_inputs no existe"}
            
        molecule_dirs = [d for d in self.base_inp_dir.iterdir() if d.is_dir()]
        
        if not molecule_dirs:
            return {"success": False, "error": "No hay moléculas para procesar"}
            
        all_results = {
            "success": True,
            "molecules": {},
            "summary": {
                "total": len(molecule_dirs),
                "completed": 0,
                "failed": 0
            }
        }
        
        print(f"🎯 Procesando {len(molecule_dirs)} moléculas...")
        
        for mol_dir in molecule_dirs:
            mol_name = mol_dir.name
            print(f"\n{'='*50}")
            print(f"🧬 Procesando molécula: {mol_name}")
            print(f"{'='*50}")
            
            results = self.process_molecule(mol_name)
            all_results["molecules"][mol_name] = results
            
            if results["success"]:
                all_results["summary"]["completed"] += 1
            else:
                all_results["summary"]["failed"] += 1
                all_results["success"] = False
                
        print(f"\n{'='*50}")
        print(f"📊 RESUMEN FINAL:")
        print(f"✅ Completadas: {all_results['summary']['completed']}")
        print(f"❌ Fallidas: {all_results['summary']['failed']}")
        print(f"{'='*50}")
        
        return all_results
        
    def get_available_molecules(self) -> List[str]:
        """Retorna lista de moléculas disponibles para procesar."""
        if not self.base_inp_dir.exists():
            return []
            
        return [d.name for d in self.base_inp_dir.iterdir() if d.is_dir()]
        
    def get_output_files(self, molecule_name: str) -> dict:
        """Retorna rutas de archivos de salida para una molécula."""
        molecule_out_dir = self.base_out_dir / molecule_name
        
        if not molecule_out_dir.exists():
            return {}
            
        output_files = {}
        calc_types = ["opt", "ir-raman", "nmr"]
        
        for calc_type in calc_types:
            out_file = molecule_out_dir / f"{molecule_name}-{calc_type}.out"
            if out_file.exists():
                output_files[calc_type] = str(out_file)
                
        return output_files
    
    def get_spectra_data(self, molecule_name: str) -> Dict[str, Any]:
        """Retorna los datos de espectros IR y Raman para una molécula."""
        molecule_out_dir = self.base_out_dir / molecule_name
        ir_raman_file = molecule_out_dir / f"{molecule_name}-ir-raman.out"
        
        if ir_raman_file.exists():
            return self.spectrum_extractor.extract_spectra_from_file(ir_raman_file)
        return {'ir_spectrum': [], 'raman_spectrum': []}


def main():
    """Función principal para usar desde línea de comandos."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generar archivos .out de ORCA y extraer espectros")
    parser.add_argument("--molecule", "-m", help="Procesar molécula específica")
    parser.add_argument("--all", "-a", action="store_true", help="Procesar todas las moléculas")
    parser.add_argument("--list", "-l", action="store_true", help="Listar moléculas disponibles")
    parser.add_argument("--spectra", "-s", help="Mostrar espectros de molécula específica")
    
    args = parser.parse_args()
    
    generator = OrcaOutputGenerator()
    
    if args.list:
        molecules = generator.get_available_molecules()
        if molecules:
            print("Moléculas disponibles:")
            for mol in molecules:
                print(f"  - {mol}")
        else:
            print("No hay moléculas disponibles")
            
    elif args.spectra:
        spectra = generator.get_spectra_data(args.spectra)
        print(f"\n📊 Espectros para {args.spectra}:")
        print("\nIR Spectrum:")
        for peak in spectra['ir_spectrum']:
            print(f"  Mode {peak['mode']}: {peak['frequency']:.2f} cm⁻¹, "
                  f"Intensity: {peak['intensity']:.2f} km/mol")
        
        print("\nRaman Spectrum:")
        for peak in spectra['raman_spectrum']:
            print(f"  Mode {peak['mode']}: {peak['frequency']:.2f} cm⁻¹, "
                  f"Activity: {peak['activity']:.6f}")
            
    elif args.molecule:
        results = generator.process_molecule(args.molecule)
        print(f"\nResultados: {results}")
        
    elif args.all:
        results = generator.process_all_molecules()
        print(f"\nResultados finales: {results['summary']}")
        
    else:
        parser.print_help()


if __name__ == "__main__":
    main()