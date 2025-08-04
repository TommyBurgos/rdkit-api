# Imagen base con Python y RDKit preinstalado
FROM intheskies/rdkit:latest

# Establece el directorio de trabajo
WORKDIR /app

# Copia el código fuente
COPY . /app

# Instala dependencias
RUN pip install --upgrade pip && \
    pip install flask gunicorn

# Expone el puerto en el que Flask corre
EXPOSE 8080

# Comando para lanzar la app con gunicorn
CMD ["gunicorn", "-w", "2", "-b", "0.0.0.0:8080", "app:app"]
