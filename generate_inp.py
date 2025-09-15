import os
import sys

class OrcaInputGenerator:
    """
    Generador de archivos .inp para ORCA a partir de un archivo .xyz.
    Crea 3 tipos de inputs: optimización, IR+Raman y RMN.
    """

    def __init__(self, xyz_path: str, base_dir: str = "orca_inputs"):
        self.xyz_path = xyz_path
        self.base_dir = base_dir
        self.molecule_name = os.path.splitext(os.path.basename(xyz_path))[0]

        # Carpeta específica para la molécula
        self.output_dir = os.path.join(self.base_dir, self.molecule_name)
        os.makedirs(self.output_dir, exist_ok=True)

        # Leer coordenadas XYZ
        self.xyz_content = self._load_xyz()

    def _load_xyz(self):
        """Carga el contenido del archivo XYZ (ignora cabecera)."""
        with open(self.xyz_path, "r") as f:
            return f.readlines()[2:]

    def _write_inp(self, filename: str, header: str):
        """Genera un archivo .inp con encabezado y coordenadas."""
        file_path = os.path.join(self.output_dir, filename)
        with open(file_path, "w") as f:
            f.write(header)
            f.write("\n* xyz 0 1\n")  # neutro, singlete
            f.writelines(self.xyz_content)
            f.write("*\n")
        print(f"✅ Archivo generado: {file_path}")

    def generate_opt_input(self):
        """Genera archivo para optimización geométrica."""
        header = "! B3LYP def2-SVP Opt Freq TightSCF\n"
        filename = f"{self.molecule_name}-opt.inp"
        self._write_inp(filename, header)

    def generate_ir_raman_input(self):
        """Genera archivo para IR + Raman."""
        header = "! B3LYP 6-31+G(d,p) NumFreq\n%elprop Polar 1 end\n"
        filename = f"{self.molecule_name}-ir-raman.inp"
        self._write_inp(filename, header)

    def generate_nmr_input(self):
        """Genera archivo para RMN."""
        header = "! B3LYP 6-311+G(2d,p) AutoAux DEFGRID2 TIGHTSCF D3BJ\n"
        filename = f"{self.molecule_name}-nmr.inp"
        self._write_inp(filename, header)

    def generate_all(self):
        """Genera los 3 inputs de una sola vez."""
        self.generate_opt_input()
        self.generate_ir_raman_input()
        self.generate_nmr_input()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python generate_inp.py nombre_molecula.xyz")
        sys.exit(1)
    xyz_file = sys.argv[1]
    # Buscar en la carpeta moleculas_xyz si no se proporciona ruta
    if not os.path.isfile(xyz_file):
        posible_path = os.path.join("moleculas_xyz", xyz_file)
        if os.path.isfile(posible_path):
            xyz_file = posible_path
        else:
            print(f"No se encontró el archivo {xyz_file} ni en moleculas_xyz/")
            sys.exit(1)
    generator = OrcaInputGenerator(xyz_file)
    generator.generate_all()
