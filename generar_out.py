import os
import sys
import subprocess
import shutil
from pathlib import Path
from typing import List, Optional
import time

class OrcaOutputGenerator:
    """
    Generador de archivos .out de ORCA a partir de archivos .inp existentes.
    Ejecuta cálculos ORCA secuencialmente: Optimización, IR/Raman, RMN.
    """

    def __init__(self, base_inp_dir: str = "orca_inputs", base_out_dir: str = "orca_outputs"):
        self.base_inp_dir = Path(base_inp_dir)
        self.base_out_dir = Path(base_out_dir)
        
        # Verificar que ORCA esté disponible
        self.orca_path = self._find_orca()
        
        # Crear carpeta de salida
        self.base_out_dir.mkdir(exist_ok=True)
        
    def _find_orca(self) -> Optional[str]:
        """Busca el ejecutable de ORCA en el sistema."""
        # Posibles ubicaciones de ORCA
        possible_paths = [
            "/opt/orca/orca",
            "/usr/local/bin/orca", 
            "/usr/bin/orca",
            "orca"  # Si está en PATH
        ]
        
        for path in possible_paths:
            if shutil.which(path):
                print(f"✅ ORCA encontrado en: {path}")
                return path
                
        print("⚠️  ORCA no encontrado. Asegúrate de que esté instalado y en PATH.")
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
            
            # Ejecutar ORCA
            cmd = [self.orca_path, str(inp_copy)]
            
            with open(out_file, "w") as f:
                process = subprocess.run(
                    cmd,
                    cwd=out_dir,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600  # 1 hora máximo por cálculo
                )
            
            if process.returncode == 0:
                print(f"✅ Cálculo completado: {out_file}")
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
            "calculations": {}
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


def main():
    """Función principal para usar desde línea de comandos."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generar archivos .out de ORCA")
    parser.add_argument("--molecule", "-m", help="Procesar molécula específica")
    parser.add_argument("--all", "-a", action="store_true", help="Procesar todas las moléculas")
    parser.add_argument("--list", "-l", action="store_true", help="Listar moléculas disponibles")
    
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
