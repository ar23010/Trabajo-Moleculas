import os
import subprocess
import shutil
import re
from pathlib import Path
from typing import List, Optional, Dict, Any
import time
import datetime

# Ruta del ejecutable de ORCA (personalizable por el usuario)
ORCA_BIN = "/home/jonathan/Trabajos_de_Orca/orca-6.1.0-f.0_linux_x86-64/bin/orca"

class OrcaSpectrumExtractor:
    """Extrae IR y Raman por separado"""

    @staticmethod
    def extract_ir(content: str) -> list[dict]:
        """Sólo IR spectrum"""
        ir_spectrum = []
        ir_pattern = r"IR SPECTRUM\s+-+\s+Mode.*?\n-+.*?\n(.*?)(?:\n\n|\Z)"
        m = re.search(ir_pattern, content, re.DOTALL)
        if m:
            lines = m.group(1).strip().splitlines()
            for line in lines:
                if line.strip() and ':' in line:
                    parts = line.split()
                    if len(parts) >= 8:
                        mode = parts[0].rstrip(':')
                        freq = parts[1]
                        eps = parts[2]
                        intensity = parts[3]
                        t_sq = parts[4]
                        tx, ty, tz = "0","0","0"
                        v = re.search(r'\((.*?)\)', line)
                        if v:
                            vec = v.group(1).split()
                            if len(vec) >= 3:
                                tx, ty, tz = vec[:3]
                        ir_spectrum.append({
                            "mode": mode, "frequency": freq, "epsilon": eps,
                            "intensity": intensity, "t_squared": t_sq,
                            "tx": tx, "ty": ty, "tz": tz})
        return ir_spectrum

    @staticmethod
    def extract_raman(content: str) -> list[dict]:
        """Sólo Raman spectrum"""
        raman_spectrum = []
        raman_pattern = r"RAMAN SPECTRUM\s+-+\s+Mode.*?\n-+.*?\n(.*?)(?:\n\n|\Z)"
        m = re.search(raman_pattern, content, re.DOTALL)
        if m:
            lines = m.group(1).strip().splitlines()
            for line in lines:
                if line.strip() and ':' in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        mode = parts[0].rstrip(':')
                        freq = parts[1]
                        activity = parts[2]
                        depol = parts[3]
                        raman_spectrum.append({
                            "mode": mode, "frequency": freq,
                            "activity": activity, "depolarization": depol})
        return raman_spectrum

    @staticmethod
    def extract_scf_energy(content: str) -> Optional[str]:
        """Extrae la última tabla completa de energía SCF total"""
        # Patrón que captura desde "TOTAL SCF ENERGY" hasta "Virial Ratio"
        scf_pattern = r"TOTAL SCF ENERGY\s+-+\s+(.*?Virial Ratio\s+:\s+[\d\.]+\s*)"
        matches = re.findall(scf_pattern, content, re.DOTALL)
        
        # Si no encuentra con el patrón completo, intenta con uno más flexible
        if not matches:
            scf_pattern = r"TOTAL SCF ENERGY\s+-+\s+(.*?)(?=\n\n|\n\s*\n|$)"
            matches = re.findall(scf_pattern, content, re.DOTALL)
        
        # Tomar la última ocurrencia (la más reciente)
        if matches:
            last_match = matches[-1]
            # Limpiar espacios en blanco extra al final
            last_match = last_match.strip()
            
            # Reconstruir la tabla manteniendo el formato original
            table = "TOTAL SCF ENERGY\n"
            table += "-" * 16 + "\n\n"
            table += last_match
            return table
        
        return None
    
    @staticmethod
    def extract_orbital_energies(content: str) -> Optional[str]:
        """Extrae la tabla de energías de orbitales"""
        # Patrón que captura la tabla completa de energías de orbitales
        orbital_pattern = r"ORBITAL ENERGIES\s+-+\s+(.*?)(?=\n\n|\n\s*\n|\Z)"
        matches = re.findall(orbital_pattern, content, re.DOTALL)

        # Tomar la última ocurrencia (la más reciente)
        if matches:
            last_match = matches[-1]
            # Limpiar y formatear la tabla
            lines = last_match.strip().splitlines()

            # Reconstruir la tabla manteniendo el formato original
            table = "ORBITAL ENERGIES\n"
            table += "-" * 17 + "\n\n"
            table += "\n".join(lines)
            return table

        return None
    
    @staticmethod
    def extract_population_analysis(content: str) -> dict:
        """Extrae los datos de población de Mulliken y Loewdin"""
        result = {
            "mulliken": {"atomic_charges": {}, "orbital_charges": {}},
            "loewdin": {"atomic_charges": {}, "orbital_charges": {}}
        }

        # Extraer MULLIKEN
        mulliken_pattern = r"MULLIKEN ATOMIC CHARGES\s+-+\s+(.*?)(?=Sum of atomic charges)"
        mulliken_match = re.search(mulliken_pattern, content, re.DOTALL)

        if mulliken_match:
            lines = mulliken_match.group(1).strip().splitlines()
            for line in lines:
                if ':' in line:
                    parts = line.split(':')
                    atom = parts[0].strip()
                    charge = float(parts[1].strip())
                    result["mulliken"]["atomic_charges"][atom] = charge

        # Extraer orbitales reducidos de MULLIKEN
        mulliken_orbital_pattern = r"MULLIKEN REDUCED ORBITAL CHARGES\s+-+\s+(.*?)(?=\n\n|\Z)"
        mulliken_orbital_match = re.search(mulliken_orbital_pattern, content, re.DOTALL)

        if mulliken_orbital_match:
            orbital_content = mulliken_orbital_match.group(1)
            atom_blocks = re.split(r'\n(?=\s*\d+\s+\w\s+)', orbital_content)

            for block in atom_blocks:
                if not block.strip():
                    continue

                lines = block.strip().splitlines()
                if not lines:
                    continue

                # Primera línea contiene el átomo
                first_line = lines[0].split(':')
                if len(first_line) >= 1:
                    atom = first_line[0].strip()
                    result["mulliken"]["orbital_charges"][atom] = {
                        "orbitals": {},
                        "s_total": 0.0,
                        "p_total": 0.0,
                        "d_total": 0.0
                    }

                    # Procesar líneas restantes
                    for line in lines[1:]:
                        line = line.strip()
                        if not line:
                            continue

                        # Buscar totales s, p, d
                        if 's :' in line:
                            parts = line.split('s :')
                            if len(parts) > 1:
                                result["mulliken"]["orbital_charges"][atom]["s_total"] = float(parts[1].strip())
                        elif 'p :' in line:
                            parts = line.split('p :')
                            if len(parts) > 1:
                                result["mulliken"]["orbital_charges"][atom]["p_total"] = float(parts[1].strip())
                        elif 'd :' in line:
                            parts = line.split('d :')
                            if len(parts) > 1:
                                result["mulliken"]["orbital_charges"][atom]["d_total"] = float(parts[1].strip())
                        else:
                            # Orbitales individuales
                            if ':' in line:
                                parts = line.split(':')
                                if len(parts) >= 2:
                                    orbital = parts[0].strip()
                                    value = float(parts[1].strip())
                                    result["mulliken"]["orbital_charges"][atom]["orbitals"][orbital] = value

        # Extraer LOEWDIN (misma estructura que Mulliken)
        loewdin_pattern = r"LOEWDIN ATOMIC CHARGES\s+-+\s+(.*?)(?=\n\n|\Z)"
        loewdin_match = re.search(loewdin_pattern, content, re.DOTALL)

        if loewdin_match:
            lines = loewdin_match.group(1).strip().splitlines()
            for line in lines:
                if ':' in line:
                    parts = line.split(':')
                    atom = parts[0].strip()
                    charge = float(parts[1].strip())
                    result["loewdin"]["atomic_charges"][atom] = charge

        # Extraer orbitales reducidos de LOEWDIN
        loewdin_orbital_pattern = r"LOEWDIN REDUCED ORBITAL CHARGES\s+-+\s+(.*?)(?=\n\n|\Z)"
        loewdin_orbital_match = re.search(loewdin_orbital_pattern, content, re.DOTALL)

        if loewdin_orbital_match:
            orbital_content = loewdin_orbital_match.group(1)
            atom_blocks = re.split(r'\n(?=\s*\d+\s+\w\s+)', orbital_content)

            for block in atom_blocks:
                if not block.strip():
                    continue

                lines = block.strip().splitlines()
                if not lines:
                    continue

                # Primera línea contiene el átomo
                first_line = lines[0].split(':')
                if len(first_line) >= 1:
                    atom = first_line[0].strip()
                    result["loewdin"]["orbital_charges"][atom] = {
                        "orbitals": {},
                        "s_total": 0.0,
                        "p_total": 0.0,
                        "d_total": 0.0
                    }

                    # Procesar líneas restantes
                    for line in lines[1:]:
                        line = line.strip()
                        if not line:
                            continue

                        # Buscar totales s, p, d
                        if 's :' in line:
                            parts = line.split('s :')
                            if len(parts) > 1:
                                result["loewdin"]["orbital_charges"][atom]["s_total"] = float(parts[1].strip())
                        elif 'p :' in line:
                            parts = line.split('p :')
                            if len(parts) > 1:
                                result["loewdin"]["orbital_charges"][atom]["p_total"] = float(parts[1].strip())
                        elif 'd :' in line:
                            parts = line.split('d :')
                            if len(parts) > 1:
                                result["loewdin"]["orbital_charges"][atom]["d_total"] = float(parts[1].strip())
                        else:
                            # Orbitales individuales
                            if ':' in line:
                                parts = line.split(':')
                                if len(parts) >= 2:
                                    orbital = parts[0].strip()
                                    value = float(parts[1].strip())
                                    result["loewdin"]["orbital_charges"][atom]["orbitals"][orbital] = value

        return result

    def extract_spectra_from_file(self, file_path: Path) -> dict:
        """Extrae IR, Raman, SCF, orbital energies y población de un archivo .out completo"""
        if not Path(file_path).exists():
            return {
                'ir_spectrum': [], 
                'raman_spectrum': [], 
                'scf_energy': None, 
                'orbital_energies': None,
                'population_data': None
            }
        content = Path(file_path).read_text(encoding="utf-8", errors="ignore")
        return {
            'ir_spectrum': self.extract_ir(content),
            'raman_spectrum': self.extract_raman(content),
            'scf_energy': self.extract_scf_energy(content),
            'orbital_energies': self.extract_orbital_energies(content),
            'population_data': self.extract_population_analysis(content)
        }


class OrcaOutputGenerator:
    """
    Generador de archivos .out de ORCA a partir de archivos .inp existentes.
    """

    def __init__(self, base_inp_dir: str = "orca_inputs", base_out_dir: str = "orca_outputs"):
        self.base_inp_dir = Path(base_inp_dir)
        self.base_out_dir = Path(base_out_dir)
        self.spectrum_extractor = OrcaSpectrumExtractor()
        self.orca_path = self._find_orca()
        self.base_out_dir.mkdir(exist_ok=True)

    def _find_orca(self) -> Optional[str]:
        """Busca el ejecutable de ORCA en el sistema o usa ruta personalizada."""
        if ORCA_BIN:
            if os.path.isfile(ORCA_BIN) and os.access(ORCA_BIN, os.X_OK):
                print(f"✅ ORCA personalizado encontrado en: {ORCA_BIN}")
                return ORCA_BIN
            else:
                print(f"❌ Ruta personalizada de ORCA no válida: {ORCA_BIN}")
                return None
        possible_paths = ["/opt/orca/orca", "/usr/local/bin/orca", "/usr/bin/orca", "orca"]
        for path in possible_paths:
            found = shutil.which(path)
            if found:
                print(f"✅ ORCA encontrado en: {found}")
                return found
        print("⚠️  ORCA no encontrado. Asegúrate de que esté instalado y en PATH o define ORCA_BIN.")
        return None

    
    def _run_orca_calculation(self, inp_file: Path, out_dir: Path) -> dict:
        """Ejecuta un cálculo ORCA individual y devuelve diccionario con info."""
        result = {"success": False, "time": 0.0}
        start = time.time()

        if not self.orca_path:
            print("❌ ORCA no está disponible")
            return result

        try:
            out_dir.mkdir(exist_ok=True)
            inp_name = inp_file.stem
            out_file = out_dir / f"{inp_name}.out"
            inp_copy = out_dir / inp_file.name
            shutil.copy2(inp_file, inp_copy)

            cmd = [self.orca_path, inp_copy.name]
            print(f"📍 Ejecutando comando: {' '.join(cmd)} en directorio: {out_dir}")

            with open(out_file, "w") as f:
                process = subprocess.run(
                    cmd,
                    cwd=str(out_dir),
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600
                )

            result["success"] = process.returncode == 0
            end = time.time()
            result["time"] = end - start

            if result["success"]:
                content = out_file.read_text(encoding="utf-8", errors="ignore")
                ir_spectrum = self.spectrum_extractor.extract_ir(content)
                raman_spectrum = self.spectrum_extractor.extract_raman(content)
                scf_energy = self.spectrum_extractor.extract_scf_energy(content)
                orbital_energies = self.spectrum_extractor.extract_orbital_energies(content)
                population_data = self.spectrum_extractor.extract_population_analysis(content)

                # Guardar individualmente
                ir_path = self._save_ir(ir_spectrum, out_dir, inp_name)
                raman_path = self._save_raman(raman_spectrum, out_dir, inp_name)
                scf_path = self._save_scf_energy(scf_energy, out_dir, inp_name)
                orbital_path = self._save_orbital_energies(orbital_energies, out_dir, inp_name)
                population_path = self._save_population_analysis(population_data, out_dir, inp_name)

                # ✅ Combinar automáticamente después de guardar los individuales
                self.combine_existing_spectra_files()

            else:
                print(f"❌ Error en cálculo: {inp_file.name}")
                print(f"Error: {process.stderr}")

            return result

        except subprocess.TimeoutExpired:
            print(f"⏰ Timeout en cálculo: {inp_file.name}")
            return result
        except Exception as e:
            print(f"❌ Error inesperado: {e}")
            return result


    def _save_ir(self, ir_spectrum, out_dir, base_name):
        if not ir_spectrum:
            return
        # Solo la carpeta moleculas (sin subcarpeta de molécula)
        modelos_dir = Path.cwd() / "modelos"
        modelos_dir.mkdir(parents=True, exist_ok=True)

        ir_file = modelos_dir / f"FINAL_ir_spectrum.txt"
        with open(ir_file, 'w') as f:
            f.write("IR SPECTRUM\n" + "-"*11 + "\n\n")
            f.write(" Mode   freq       eps      Int      T**2         TX        TY        TZ\n")
            f.write("       cm**-1   L/(mol*cm) km/mol    a.u.\n")
            f.write("-"*75 + "\n")
            for p in ir_spectrum:
                f.write(f"{p['mode']:>4}: {p['frequency']:>8}   {p['epsilon']:>8}   "
                        f"{p['intensity']:>6}   {p['t_squared']:>8}  "
                        f"( {p['tx']:>8}  {p['ty']:>8}  {p['tz']:>8})\n")
        print(f"📊 IR guardado en: {ir_file}")

    def _save_raman(self, raman_spectrum, out_dir, base_name):
        if not raman_spectrum:
            return
        modelos_dir = Path.cwd() / "modelos"
        modelos_dir.mkdir(parents=True, exist_ok=True)

        raman_file = modelos_dir / f"FINAL_raman_spectrum.txt"
        with open(raman_file, 'w') as f:
            f.write("RAMAN SPECTRUM\n" + "-"*15 + "\n\n")
            f.write(" Mode    freq (cm**-1)   Activity   Depolarization\n")
            f.write("-"*55 + "\n")
            for p in raman_spectrum:
                f.write(f"{p['mode']:>5}: {p['frequency']:>11}   "
                        f"{p['activity']:>9}   {p['depolarization']:>12}\n")
        print(f"📊 Raman guardado en: {raman_file}")

    def _save_scf_energy(self, scf_energy, out_dir, base_name):
        """Guarda la tabla de energía SCF en un archivo separado"""
        if not scf_energy:
            return
        modelos_dir = Path.cwd() / "modelos"
        modelos_dir.mkdir(parents=True, exist_ok=True)

        scf_file = modelos_dir / f"FINAL_scf_energy.txt"
        with open(scf_file, 'w') as f:
            f.write(scf_energy)
        print(f"📊 Energía SCF guardada en: {scf_file}")

    def _save_orbital_energies(self, orbital_energies, out_dir, base_name):
        """Guarda la tabla de energías orbitales en un archivo separado"""
        if not orbital_energies:
            return
        modelos_dir = Path.cwd() / "modelos"
        modelos_dir.mkdir(parents=True, exist_ok=True)

        orbital_file = modelos_dir / f"FINAL_orbital_energies.txt"
        with open(orbital_file, 'w') as f:
            f.write(orbital_energies)
        print(f"📊 Energías orbitales guardadas en: {orbital_file}")

    def _save_population_analysis(self, population_data, out_dir, base_name):
        """Guarda el análisis de población en un archivo separado"""
        if not population_data or (not population_data["mulliken"]["atomic_charges"] and not population_data["loewdin"]["atomic_charges"]):
            return

        modelos_dir = Path.cwd() / "modelos"
        modelos_dir.mkdir(parents=True, exist_ok=True)

        pop_file = modelos_dir / f"FINAL_population_analysis.txt"

        with open(pop_file, 'w') as f:
            # Escribir datos de Mulliken
            f.write("MULLIKEN POPULATION ANALYSIS\n")
            f.write("=" * 30 + "\n\n")
            f.write("Atomic Charges:\n")
            f.write("-" * 15 + "\n")
            for atom, charge in population_data["mulliken"]["atomic_charges"].items():
                f.write(f"{atom}: {charge:10.6f}\n")

            f.write("\nReduced Orbital Charges:\n")
            f.write("-" * 25 + "\n")
            for atom, data in population_data["mulliken"]["orbital_charges"].items():
                f.write(f"\n{atom}:\n")
                f.write(f"  s total: {data['s_total']:10.6f}\n")
                f.write(f"  p total: {data['p_total']:10.6f}\n")
                if data['d_total'] != 0.0:
                    f.write(f"  d total: {data['d_total']:10.6f}\n")
                for orbital, value in data["orbitals"].items():
                    f.write(f"  {orbital}: {value:10.6f}\n")

            # Escribir datos de Loewdin
            f.write("\n\nLOEWDIN POPULATION ANALYSIS\n")
            f.write("=" * 30 + "\n\n")
            f.write("Atomic Charges:\n")
            f.write("-" * 15 + "\n")
            for atom, charge in population_data["loewdin"]["atomic_charges"].items():
                f.write(f"{atom}: {charge:10.6f}\n")

            f.write("\nReduced Orbital Charges:\n")
            f.write("-" * 25 + "\n")
            for atom, data in population_data["loewdin"]["orbital_charges"].items():
                f.write(f"\n{atom}:\n")
                f.write(f"  s total: {data['s_total']:10.6f}\n")
                f.write(f"  p total: {data['p_total']:10.6f}\n")
                if data['d_total'] != 0.0:
                    f.write(f"  d total: {data['d_total']:10.6f}\n")
                for orbital, value in data["orbitals"].items():
                    f.write(f"  {orbital}: {value:10.6f}\n")

        print(f"📊 Análisis de población guardado en: {pop_file}")    

    def combine_existing_spectra_files(self):
        """
        Combina los archivos FINAL_ir_spectrum.txt y FINAL_raman_spectrum.txt
        en un solo archivo FINAL_combined_spectra.txt automáticamente
        """
        modelos_dir = Path.cwd() / "modelos"
        
        ir_file = modelos_dir / "FINAL_ir_spectrum.txt"
        raman_file = modelos_dir / "FINAL_raman_spectrum.txt"
        combined_file = modelos_dir / "FINAL_combined_spectra.txt"
        
        # Verificar que existan los archivos
        if not ir_file.exists() and not raman_file.exists():
            print("⚠️  No se encontraron archivos de espectros para combinar")
            return None
        
        try:
            with open(combined_file, 'w', encoding='utf-8') as outfile:
                # Combinar archivo IR si existe
                if ir_file.exists():
                    print(f"📥 Leyendo archivo IR: {ir_file}")
                    with open(ir_file, 'r', encoding='utf-8') as infile:
                        outfile.write(infile.read())
                        outfile.write("\n\n")  # Separador entre espectros
                
                # Combinar archivo Raman si existe
                if raman_file.exists():
                    print(f"📥 Leyendo archivo Raman: {raman_file}")
                    with open(raman_file, 'r', encoding='utf-8') as infile:
                        outfile.write(infile.read())
            
            print(f"✅ Archivos combinados exitosamente en: {combined_file}")
            return combined_file
            
        except Exception as e:
            print(f"❌ Error al combinar archivos: {e}")
            return None

    def process_molecule(self, molecule_name: str) -> dict:
        """
        Procesa todos los .inp de una molécula y devuelve:
        {
          "success": bool,
          "calculations": {
             "opt": {"success": bool, "time": float},
             "ir-raman": {"success": bool, "time": float},
             ...
          }
        }
        """
        inp_dir = self.base_inp_dir / molecule_name
        out_dir = self.base_out_dir / molecule_name
    
        if not inp_dir.exists():
            return {"success": False, "error": f"No existe carpeta {inp_dir}"}
    
        inp_files = list(inp_dir.glob("*.inp"))
        if not inp_files:
            return {"success": False, "error": f"No hay .inp en {inp_dir}"}
    
        out_dir.mkdir(exist_ok=True)
    
        calculations = {}
        for inp_file in inp_files:
            # tomar la parte después del último guion como tipo
            calc_type = inp_file.stem.split("-")[-1]  # "opt" o "ir-raman"
            result = self._run_orca_calculation(inp_file, out_dir)
            calculations[calc_type] = result
    
        success = all(c["success"] for c in calculations.values())
        return {"success": success, "calculations": calculations}


    def process_all_molecules(self) -> dict:
        """Procesa todos los cálculos ORCA para todas las moléculas disponibles."""
        if not self.base_inp_dir.exists():
            return {"success": False, "error": "Carpeta orca_inputs no existe"}
        molecule_dirs = [d for d in self.base_inp_dir.iterdir() if d.is_dir()]
        if not molecule_dirs:
            return {"success": False, "error": "No hay moléculas para procesar"}
        all_results = {"success": True, "molecules": {}, "summary": {"total": len(molecule_dirs), "completed": 0, "failed": 0}}
        print(f"🎯 Procesando {len(molecule_dirs)} moléculas...")
        for mol_dir in molecule_dirs:
            mol_name = mol_dir.name
            print(f"\n{'='*50}\n🧬 Procesando molécula: {mol_name}\n{'='*50}")
            results = self.process_molecule(mol_name)
            all_results["molecules"][mol_name] = results
            if results["success"]:
                all_results["summary"]["completed"] += 1
            else:
                all_results["summary"]["failed"] += 1
                all_results["success"] = False
        print(f"\n{'='*50}\n📊 RESUMEN FINAL:")
        print(f"✅ Completadas: {all_results['summary']['completed']}")
        print(f"❌ Fallidas: {all_results['summary']['failed']}\n{'='*50}")
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
        return {'ir_spectrum': [], 'raman_spectrum': [], 'scf_energy': None}
    
    def get_population_data(self, molecule_name: str) -> Dict[str, Any]:
        """Retorna los datos de análisis de población para una molécula."""
        molecule_out_dir = self.base_out_dir / molecule_name
        ir_raman_file = molecule_out_dir / f"{molecule_name}-ir-raman.out"
        if ir_raman_file.exists():
            content = ir_raman_file.read_text(encoding="utf-8", errors="ignore")
            return self.spectrum_extractor.extract_population_analysis(content)
        return {"mulliken": {"atomic_charges": {}, "orbital_charges": {}}, "loewdin": {"atomic_charges": {}, "orbital_charges": {}}}


def main():
    """Función principal para usar desde línea de comandos."""
    import argparse
    parser = argparse.ArgumentParser(description="Generar archivos .out de ORCA y extraer espectros")
    parser.add_argument("--molecule", "-m", help="Procesar molécula específica")
    parser.add_argument("--all", "-a", action="store_true", help="Procesar todas las moléculas")
    parser.add_argument("--list", "-l", action="store_true", help="Listar moléculas disponibles")
    parser.add_argument("--spectra", "-s", help="Mostrar espectros de molécula específica")
    parser.add_argument("--scf", help="Mostrar energía SCF de molécula específica")
    parser.add_argument("--orbitals", help="Mostrar energías orbitales de molécula específica")
    parser.add_argument("--population", "-p", help="Mostrar análisis de población de molécula específica")
    args = parser.parse_args()
    generator = OrcaOutputGenerator()
    
    if args.list:
        molecules = generator.get_available_molecules()
        if molecules:
            print("Moléculas disponibles:")
            for mol in molecules: print(f"  - {mol}")
        else:
            print("No hay moléculas disponibles")
            
    elif args.spectra:
        spectra = generator.get_spectra_data(args.spectra)
        print(f"\n📊 ESPECTROS FINALES para {args.spectra}:")
        print(f"\nIR Spectrum ({len(spectra['ir_spectrum'])} modos):")
        print("Mode   freq       eps      Int      T**2         TX        TY        TZ")
        print("----------------------------------------------------------------------------")
        for peak in spectra['ir_spectrum']:
            print(f"{peak['mode']:>4}: {peak['frequency']:>8}   {peak['epsilon']:>8}   "
                  f"{peak['intensity']:>6}   {peak['t_squared']:>8}  "
                  f"( {peak['tx']:>8}  {peak['ty']:>8}  {peak['tz']:>8})")
        print(f"\nRaman Spectrum ({len(spectra['raman_spectrum'])} modos):")
        print("Mode    freq (cm**-1)   Activity   Depolarization")
        print("-------------------------------------------------------------------")
        for peak in spectra['raman_spectrum']:
            print(f"{peak['mode']:>5}: {peak['frequency']:>11}   "
                  f"{peak['activity']:>9}   {peak['depolarization']:>12}")
                  
    elif args.scf:
        spectra = generator.get_spectra_data(args.scf)
        if spectra['scf_energy']:
            print(f"\n📊 ENERGÍA SCF para {args.scf}:")
            print(spectra['scf_energy'])
        else:
            print(f"No se encontró información de energía SCF para {args.scf}")
            
    elif args.orbitals:
        spectra = generator.get_spectra_data(args.orbitals)
        if spectra['orbital_energies']:
            print(f"\n📊 ENERGÍAS ORBITALES para {args.orbitals}:")
            print(spectra['orbital_energies'])
        else:
            print(f"No se encontró información de energías orbitales para {args.orbitals}")
            
    elif args.population:
        # Necesitamos una nueva función para obtener datos de población
        population_data = generator.get_population_data(args.population)
        if population_data and (population_data["mulliken"]["atomic_charges"] or population_data["loewdin"]["atomic_charges"]):
            print(f"\n📊 ANÁLISIS DE POBLACIÓN para {args.population}:")
            
            print("\nMULLIKEN POPULATION ANALYSIS")
            print("=" * 30)
            print("\nAtomic Charges:")
            print("-" * 15)
            for atom, charge in population_data["mulliken"]["atomic_charges"].items():
                print(f"{atom}: {charge:10.6f}")
            
            print("\nReduced Orbital Charges:")
            print("-" * 25)
            for atom, data in population_data["mulliken"]["orbital_charges"].items():
                print(f"\n{atom}:")
                print(f"  s total: {data['s_total']:10.6f}")
                print(f"  p total: {data['p_total']:10.6f}")
                if data['d_total'] != 0.0:
                    print(f"  d total: {data['d_total']:10.6f}")
                for orbital, value in data["orbitals"].items():
                    print(f"  {orbital}: {value:10.6f}")
            
            print("\n\nLOEWDIN POPULATION ANALYSIS")
            print("=" * 30)
            print("\nAtomic Charges:")
            print("-" * 15)
            for atom, charge in population_data["loewdin"]["atomic_charges"].items():
                print(f"{atom}: {charge:10.6f}")
            
            print("\nReduced Orbital Charges:")
            print("-" * 25)
            for atom, data in population_data["loewdin"]["orbital_charges"].items():
                print(f"\n{atom}:")
                print(f"  s total: {data['s_total']:10.6f}")
                print(f"  p total: {data['p_total']:10.6f}")
                if data['d_total'] != 0.0:
                    print(f"  d total: {data['d_total']:10.6f}")
                for orbital, value in data["orbitals"].items():
                    print(f"  {orbital}: {value:10.6f}")
        else:
            print(f"No se encontró información de análisis de población para {args.population}")
            
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