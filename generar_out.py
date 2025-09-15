import os
import subprocess
import shutil
import re
from pathlib import Path
from typing import List, Optional, Dict, Any
import time
import datetime

# Ruta del ejecutable de ORCA (personalizable por el usuario)
ORCA_BIN = "/home/fercho/orca_6_1_0_linux_x86-64_shared_openmpi418/orca"

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

    def extract_spectra_from_file(self, file_path: Path) -> dict:
        """Extrae IR y Raman de un archivo .out completo"""
        if not Path(file_path).exists():
            return {'ir_spectrum': [], 'raman_spectrum': []}
        content = Path(file_path).read_text(encoding="utf-8", errors="ignore")
        return {
            'ir_spectrum': self.extract_ir(content),
            'raman_spectrum': self.extract_raman(content)
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

                # Guardar individualmente
                ir_path = self._save_ir(ir_spectrum, out_dir, inp_name)
                raman_path = self._save_raman(raman_spectrum, out_dir, inp_name)
                print("depuracion")

                # Combinar automáticamente los archivos creados - VERIFICAR ESTA LÍNEA
                if ir_path and raman_path:
                    print("true")
                    combined_path = self.combinar_archivos(ir_path, raman_path)
                    if combined_path:
                        print(f"✅ Archivos combinados en: {combined_path}")
                    else:
                        print("❌ Error al combinar archivos")
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

    def combinar_archivos(self, archivo1, archivo2):
        """
        Combina dos archivos de texto y guarda el resultado en modelos/FINAL_spectrum.txt
        """
        try:
            # Convertir a Path si son strings
            if isinstance(archivo1, str):
                archivo1 = Path(archivo1)
            if isinstance(archivo2, str):
                archivo2 = Path(archivo2)

            # Verificar que ambos archivos existen
            if not archivo1.exists():
                print(f"❌ Archivo no encontrado: {archivo1}")
                return None
            if not archivo2.exists():
                print(f"❌ Archivo no encontrado: {archivo2}")
                return None

            # Crear directorio modelos si no existe
            modelos_dir = Path.cwd() / "modelos"
            modelos_dir.mkdir(parents=True, exist_ok=True)

            # Ruta del archivo combinado
            combined_file = modelos_dir / "FINAL_spectrum.txt"

            # Leer el contenido de ambos archivos
            with open(archivo1, 'r', encoding='utf-8') as f1:
                contenido1 = f1.read()

            with open(archivo2, 'r', encoding='utf-8') as f2:
                contenido2 = f2.read()

            # Combinar los contenidos con separación clara
            separador = "\n" + "="*80 + "\n"
            contenido_combinado = f"ARCHIVO COMBINADO - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            contenido_combinado += separador + "ESPECTRO IR:\n" + separador
            contenido_combinado += contenido1
            contenido_combinado += separador + "ESPECTRO RAMAN:\n" + separador
            contenido_combinado += contenido2

            # Escribir el resultado en el archivo de salida
            with open(combined_file, 'w', encoding='utf-8') as f_salida:
                f_salida.write(contenido_combinado)

            print(f"✅ Archivos combinados exitosamente en: {combined_file}")
            return combined_file

        except FileNotFoundError as e:
            print(f"❌ Error: Archivo no encontrado - {e}")
            return None
        except Exception as e:
            print(f"❌ Error inesperado al combinar: {e}")
            import traceback
            traceback.print_exc()
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
