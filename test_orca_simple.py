#!/usr/bin/env python3
"""Test simple para verificar ORCA desde Python"""

import subprocess
from pathlib import Path

# Configuración
orca_path = "/home/jonathan/Trabajos de Orca/orca-6.1.0-f.0_linux_x86-64/bin/orca"
work_dir = Path("/home/jonathan/Trabajo-Moleculas/orca_outputs/co2")
inp_file = "co2-opt.inp"
out_file = "test_python.out"

print(f"🔍 Directorio de trabajo: {work_dir}")
print(f"🔍 Archivo de entrada: {inp_file}")
print(f"🔍 Ejecutable ORCA: {orca_path}")

# Comando
cmd = [orca_path, inp_file]

print(f"🚀 Ejecutando: {' '.join(cmd)}")

try:
    with open(work_dir / out_file, "w") as f:
        process = subprocess.run(
            cmd,
            cwd=str(work_dir),  # Importante: directorio de trabajo
            stdout=f,
            stderr=subprocess.PIPE,
            text=True,
            timeout=300  # 5 minutos para test
        )
    
    print(f"📊 Código de salida: {process.returncode}")
    
    if process.returncode == 0:
        print("✅ ¡ORCA ejecutado exitosamente desde Python!")
    else:
        print("❌ Error en ORCA:")
        print(process.stderr)
        
except subprocess.TimeoutExpired:
    print("⏰ Timeout del proceso")
except Exception as e:
    print(f"❌ Error inesperado: {e}")
