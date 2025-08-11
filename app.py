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


from rdkit.Chem import AllChem, SDWriter
import tempfile

app = Flask(__name__)

@app.route("/modificar", methods=["POST"])
def modificar_molecula():
    try:
        data = request.get_json()
        smiles = data.get("smiles")
        atomo_idx = int(data.get("atomo_idx"))
        grupo = data.get("grupo")

        if not smiles or grupo is None:
            return jsonify({"error": "Datos insuficientes"}), 400

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        edmol = Chem.RWMol(mol)

        # Remover un hidrógeno unido al átomo objetivo
        h_idx = None
        for vecino in edmol.GetAtomWithIdx(atomo_idx).GetNeighbors():
            if vecino.GetAtomicNum() == 1:
                h_idx = vecino.GetIdx()
                break

        if h_idx is not None:
            edmol.RemoveAtom(h_idx)

        if grupo == "CH3":
            c = Chem.Atom(6)
            h1 = Chem.Atom(1)
            h2 = Chem.Atom(1)
            h3 = Chem.Atom(1)
            c_idx = edmol.AddAtom(c)
            edmol.AddBond(atomo_idx, c_idx, Chem.rdchem.BondType.SINGLE)
            edmol.AddBond(c_idx, edmol.AddAtom(h1), Chem.rdchem.BondType.SINGLE)
            edmol.AddBond(c_idx, edmol.AddAtom(h2), Chem.rdchem.BondType.SINGLE)
            edmol.AddBond(c_idx, edmol.AddAtom(h3), Chem.rdchem.BondType.SINGLE)

        mol_final = edmol.GetMol()
        Chem.SanitizeMol(mol_final)
        mol_final.UpdatePropertyCache(strict=False)
        AllChem.EmbedMolecule(mol_final)
        AllChem.UFFOptimizeMolecule(mol_final)

        # Guardar en SDF temporal
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp:
            writer = SDWriter(temp.name)
            writer.write(mol_final)
            writer.close()
            with open(temp.name, "r") as f:
                sdf_data = f.read()

        return jsonify({"sdf": sdf_data})

    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/")
def home():
    return "✅ Microservicio RDKit activo"
