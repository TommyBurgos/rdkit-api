# Imagen base oficial con RDKit y Python
FROM rdkit/rdkit:2022.09.5-py3.9

# Establece el directorio de trabajo
WORKDIR /app

# Copia el código fuente
COPY . /app

# Instala dependencias adicionales
RUN pip install --upgrade pip && \
    pip install flask gunicorn

# Expone el puerto que usará la app
EXPOSE 8080

# Comando para correr el servidor
CMD ["gunicorn", "-w", "2", "-b", "0.0.0.0:8080", "app:app"]
