FROM continuumio/miniconda3

# Crear y activar entorno con rdkit
RUN conda create -n rdkit_env -c rdkit -c conda-forge rdkit flask gunicorn python=3.9 -y

# Activar el entorno
SHELL ["conda", "run", "-n", "rdkit_env", "/bin/bash", "-c"]

# Copiar la app
WORKDIR /app
COPY . .

# Exponer puerto
EXPOSE 8080

# Ejecutar gunicorn usando el entorno conda
CMD ["conda", "run", "-n", "rdkit_env", "gunicorn", "-b", "0.0.0.0:8080", "app:app"]
