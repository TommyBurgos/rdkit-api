from flask import Flask, request, jsonify
from rdkit import Chem

app = Flask(__name__)

@app.route("/")
def inicio():
    return "¡API RDKit funcionando!"

@app.route("/validar", methods=["POST"])
def validar_molecula():
    data = request.get_json()
    smiles = data.get("smiles")
    objetivo = data.get("objetivo")

    if not smiles:
        return jsonify({"error": "SMILES no proporcionado"}), 400

    try:
        mol_estudiante = Chem.MolFromSmiles(smiles)
        if mol_estudiante is None:
            return jsonify({"valida": False, "razon": "SMILES inválido"})

        if objetivo:
            mol_objetivo = Chem.MolFromSmiles(objetivo)
            es_igual = Chem.MolToSmiles(mol_estudiante) == Chem.MolToSmiles(mol_objetivo)
        else:
            es_igual = None

        return jsonify({"valida": True, "coincide_objetivo": es_igual})
    except Exception as e:
        return jsonify({"valida": False, "razon": str(e)})

# 👇 Esto es lo que faltaba para ejecutarlo
if __name__ == "__main__":
    app.run(debug=True)
