# 🧬 Analizador de Moléculas XYZ con ORCA

Aplicación web desarrollada en Streamlit para el análisis de moléculas utilizando archivos XYZ y cálculos cuánticos con ORCA.

## 📋 Características

### ✅ Funciones Implementadas

1. **🔄 Generación automática de archivos ORCA (.inp)**
   - Optimización geométrica
   - Espectros IR/Raman 
   - RMN
   - Usa funcional B3LYP con bases 6-31+G(d,p) y 6-311+G(2d,p)

2. **⚙️ Ejecución de cálculos ORCA**
   - Ejecuta cálculos secuencialmente
   - Crea carpeta `orca_outputs` automáticamente
   - Manejo de timeouts y errores

3. **📊 Interfaz interactiva mejorada**
   - Selector de moléculas desde carpeta `moleculas_xyz`
   - Información detallada de cada molécula
   - Visualizaciones 3D y 2D
   - Análisis espectroscópico

4. **🧪 Funciones de análisis**
   - Comparación con moléculas de referencia (NH₃)
   - Función de distribución radial (RDF)
   - Espectros IR y Raman
   - Análisis de desplazamientos químicos (RMN)

## 📁 Estructura del Proyecto

```
Trabajo-Moleculas/
├── spectra_app.py          # Aplicación principal de Streamlit
├── generate_inp.py         # Generador de archivos .inp para ORCA
├── generar_out.py         # Ejecutor de cálculos ORCA
├── moleculas_xyz/         # Carpeta con archivos .xyz
│   ├── benzene.xyz
│   ├── water.xyz
│   ├── ethanol.xyz
│   └── ...
├── orca_inputs/           # Archivos .inp generados
└── orca_outputs/          # Resultados de ORCA (.out)
```

## 🚀 Instalación y Uso

### Requisitos

1. **Python 3.8+** con las siguientes librerías:
   ```bash
   pip install streamlit numpy matplotlib plotly pandas pathlib
   ```

2. **ORCA** (opcional, para cálculos cuánticos):
   - Descargar e instalar desde [ORCA Forum](https://orcaforum.kofo.mpg.de)
   - Asegurar que `orca` esté en el PATH del sistema

### Ejecución

1. **Clonar o descargar el proyecto**
2. **Navegar a la carpeta del proyecto**:
   ```bash
   cd Trabajo-Moleculas
   ```

3. **Ejecutar la aplicación**:
   ```bash
   streamlit run spectra_app.py
   ```

4. **Abrir en el navegador**: La aplicación se abrirá automáticamente en `http://localhost:8501`

## 📖 Guía de Uso

### 1. Selección de Molécula
- Al abrir la aplicación, el menú lateral mostrará todas las moléculas disponibles en `moleculas_xyz/`
- Seleccionar la molécula que se desea analizar

### 2. Tipos de Análisis Disponibles

#### 📋 **Información de la molécula**
- Muestra datos básicos: número de átomos, elementos, fórmula
- Tabla de coordenadas atómicas
- Descripción del archivo XYZ

#### 🔄 **Procesar con ORCA**
- **Paso 1**: Genera automáticamente archivos .inp para:
  - Optimización geométrica (`molecula-opt.inp`)
  - IR/Raman (`molecula-ir-raman.inp`)
  - RMN (`molecula-nmr.inp`)
- **Paso 2**: Ejecuta cálculos ORCA (si está disponible)
- **Paso 3**: Guarda resultados en `orca_outputs/`

#### 🧪 **Visualizaciones**
- **3D/2D**: Representaciones moleculares interactivas
- **Contenedor**: Simulación de múltiples moléculas
- **Espectros**: IR, Raman, RMN

#### 🔍 **Análisis Comparativo**
- Comparación con moléculas de referencia
- Análisis de diferencias geométricas
- Métricas moleculares

## 🛠️ Configuración Avanzada

### Personalizar Cálculos ORCA

En `generate_inp.py`, puedes modificar los parámetros de cálculo:

```python
# Optimización
header = "! B3LYP 6-31+G(d,p) OPT DEFGRID2 TIGHTSCF D3BJ PAL8\n"

# IR/Raman  
header = "! B3LYP 6-31+G(d,p) NumFreq\n%elprop Polar 1 end\n"

# RMN
header = "! B3LYP 6-311+G(2d,p) AutoAux DEFGRID2 TIGHTSCF D3BJ\n"
```

### Agregar Nuevas Moléculas

1. Crear archivo `.xyz` en formato estándar:
   ```
   3
   Descripción de la molécula
   O    0.000000    0.000000    0.000000
   H    0.757000    0.586000    0.000000  
   H   -0.757000    0.586000    0.000000
   ```

2. Guardar en la carpeta `moleculas_xyz/`
3. La aplicación detectará automáticamente el nuevo archivo

## 🔧 Solución de Problemas

### Errores Comunes

1. **"ORCA no encontrado"**
   - Instalar ORCA y agregarlo al PATH
   - Los archivos .inp se generan independientemente

2. **"Carpeta moleculas_xyz no existe"**
   - Crear la carpeta y agregar archivos .xyz

3. **Errores de dependencias**
   - Instalar todas las librerías requeridas: `pip install streamlit numpy matplotlib plotly pandas`

### Logs y Debugging

- Los errores de ORCA se muestran en la consola
- Los archivos .inp se pueden ejecutar manualmente si ORCA falla
- Verificar permisos de escritura en las carpetas de salida

## 📊 Formatos de Salida

### Archivos Generados

- **`.inp`**: Archivos de entrada para ORCA
- **`.out`**: Resultados de cálculos ORCA
- **`.gbw`, `.opt`, `.freq`**: Archivos auxiliares de ORCA

### Datos Exportables

- Tablas de coordenadas (CSV/Excel desde Streamlit)
- Gráficos (PNG desde matplotlib/plotly)
- Métricas moleculares (DataFrame de pandas)

## 👥 Desarrollo

### Estructura del Código

- **`spectra_app.py`**: Interfaz principal y funciones de análisis
- **`generate_inp.py`**: Clase `OrcaInputGenerator` para crear inputs
- **`generar_out.py`**: Clase `OrcaOutputGenerator` para ejecutar ORCA

### Agregar Nuevas Funciones

1. Implementar la función en `spectra_app.py`
2. Agregar opción al menú en `main()`
3. Documentar en este README

## 📄 Licencia

Este proyecto fue desarrollado para análisis científico de moléculas. Uso libre para fines académicos y de investigación.

## 🤝 Contribuir

Si encuentras bugs o tienes sugerencias:
1. Crear un issue describiendo el problema
2. Proponer mejoras o nuevas funcionalidades
3. Enviar pull requests con código limpio y documentado

---

**Desarrollado por:** Experto en Ciencias Computacionales  
**Fecha:** Septiembre 2025  
**Versión:** 2.0
