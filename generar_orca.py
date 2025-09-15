import os

class OrcaInputGenerator:
    """
    Generador de archivos .inp para ORCA a partir de un archivo .xyz.
    Crea 3 tipos de inputs: optimización, IR+Raman y RMN.
    """

    def __init__(self, xyz_path: str, output_dir: str = "orca_inputs"):
        self.xyz_path = xyz_path
        self.output_dir = output_dir
        self.molecule_name = os.path.splitext(os.path.basename(xyz_path))[0]

        # Leer coordenadas XYZ
        self.xyz_content = self._load_xyz()

        # Crear carpeta de salida
        os.makedirs(self.output_dir, exist_ok=True)

    def _load_xyz(self):
        """Carga el contenido del archivo XYZ."""
        with open(self.xyz_path, "r") as f:
            return f.readlines()[2:]  # ignoramos num átomos y comentario

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
        header = "! B3LYP 6-31+G(d,p) OPT DEFGRID2 TIGHTSCF D3BJ PAL8\n"
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


# Ejemplo de uso:
if __name__ == "__main__":
    xyz_file = "acetanilide.xyz"  # aquí va el archivo exportado desde Avogadro
    generator = OrcaInputGenerator(xyz_file)
    generator.generate_all()
