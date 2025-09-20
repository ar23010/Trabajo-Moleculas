import matplotlib.pyplot as plt
import numpy as np

# Datos de energía en Eh (Hartree)
energias_eh = {
    'Total Energy': -76.32126578130973,
    'Nuclear Repulsion': 9.10440713518709,
    'Electronic Energy': -85.42571817874006,
    'One Electron Energy': -122.95596269599451,
    'Two Electron Energy': 37.53024451725446,
    'Potential Energy': -152.15284911004716,
    'Kinetic Energy': 75.83158332873744
}

# Crear figura con 2 gráficos principales
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Análisis Energético de la Molécula', fontsize=16, fontweight='bold')

# 1. Gráfico de componentes principales (Eh)
categorias = ['Nuclear', 'One e⁻', 'Two e⁻', 'Electronic', 'Total']
valores_eh = [
    energias_eh['Nuclear Repulsion'],
    energias_eh['One Electron Energy'],
    energias_eh['Two Electron Energy'],
    energias_eh['Electronic Energy'],
    energias_eh['Total Energy']
]

colors = ['#FF6B6B', '#45B7D1', '#F9A602', '#4ECDC4', '#2E86AB']
bars = ax1.bar(categorias, valores_eh, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
ax1.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
ax1.set_ylabel('Energía (Eh)', fontweight='bold')
ax1.set_title('Componentes de Energía Molecular (Hartree)', fontweight='bold')
ax1.grid(axis='y', alpha=0.3, linestyle='--')

# Añadir valores en las barras (Eh)
for bar, valor in zip(bars, valores_eh):
    height = bar.get_height()
    va_position = 'bottom' if height > 0 else 'top'
    y_offset = 0.8 if height > 0 else -1.2
    ax1.text(bar.get_x() + bar.get_width()/2., height + y_offset,
             f'{valor:.2f} Eh', ha='center', va=va_position, fontsize=9, fontweight='bold')

# 2. Gráfico del teorema del virial
virial_categorias = ['Energía Potencial', 'Energía Cinética', 'Energía Total']
virial_valores_eh = [
    energias_eh['Potential Energy'],
    energias_eh['Kinetic Energy'],
    energias_eh['Total Energy']
]

bars_virial = ax2.bar(virial_categorias, virial_valores_eh, 
                     color=['#FF6B6B', '#4ECDC4', '#2E86AB'], alpha=0.8, 
                     edgecolor='black', linewidth=0.5)

ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
ax2.set_ylabel('Energía (Eh)', fontweight='bold')
ax2.set_title('Teorema del Virial', fontweight='bold')
ax2.grid(axis='y', alpha=0.3, linestyle='--')

# Añadir valores en las barras del virial
for bar, valor in zip(bars_virial, virial_valores_eh):
    height = bar.get_height()
    va_position = 'bottom' if height > 0 else 'top'
    y_offset = 3 if height > 0 else -4
    ax2.text(bar.get_x() + bar.get_width()/2., height + y_offset,
             f'{valor:.2f} Eh', ha='center', va=va_position, fontsize=9, fontweight='bold')

# Añadir ratio virial
virial_ratio = abs(energias_eh['Potential Energy'] / energias_eh['Kinetic Energy'])
ax2.text(0.95, 0.95, f'Ratio Virial: {virial_ratio:.4f}\n(ideal ≈ 2.0)',
         transform=ax2.transAxes, ha='right', va='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()

# Mostrar resumen numérico completo
print("="*60)
print("RESUMEN ENERGÉTICO COMPLETO - MOLÉCULA")
print("="*60)
print(f"{'COMPONENTE':<25} {'ENERGÍA (Eh)':<20} {'TIPO':<15}")
print("-"*60)
for key in energias_eh.keys():
    tipo = "Nuclear" if "Nuclear" in key else "Electrónica" if "Electronic" in key else "Mixta"
    print(f"{key:<25} {energias_eh[key]:<20.6f} {tipo:<15}")
print("="*60)
print(f"Ratio Virial (|V/T|): {virial_ratio:.6f}")
print("="*60)