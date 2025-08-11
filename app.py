from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

import json
import re
import urllib.request
import urllib.parse

app = Flask(__name__)

@app.route("/")
def home():
    return "✅ Microservicio RDKit activo"

# -------- EXISTENTES (déjalas si ya las usas) --------
@app.route("/validar", methods=["POST"])
def validar_molecula():
    data = request.get_json(silent=True) or {}
    smiles = data.get("smiles")
    objetivo = data.get("objetivo")
    if not smiles:
        return jsonify({"error": "SMILES no proporcionado"}), 400

    try:
        mol_est = Chem.MolFromSmiles(smiles)
        if mol_est is None:
            return jsonify({"valida": False, "razon": "SMILES inválido"})
        if objetivo:
            mol_obj = Chem.MolFromSmiles(objetivo)
            es_igual = Chem.MolToSmiles(mol_est, canonical=True) == Chem.MolToSmiles(mol_obj, canonical=True)
        else:
            es_igual = None
        return jsonify({"valida": True, "coincide_objetivo": es_igual})
    except Exception as e:
        return jsonify({"valida": False, "razon": str(e)}), 500

from rdkit.Chem import AllChem, SDWriter
import tempfile

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


# Si ya tienes /procesar en otro archivo, OK; si no, aquí un stub compatible:
@app.route("/procesar", methods=["POST"])
def procesar():
    """
    Espera: { "entrada": <smiles|mol>, "tipo": "smiles"|"mol" }
    Devuelve (ejemplo): {"smiles_canonico":"...", "sdf":"..."} si quieres.
    """
    data = request.get_json(silent=True) or {}
    entrada = (data.get("entrada") or "").strip()
    tipo = (data.get("tipo") or "smiles").lower()
    if not entrada:
        return jsonify({"error":"entrada vacía"}), 400
    try:
        if tipo == "mol":
            mol = Chem.MolFromMolBlock(entrada, sanitize=True)
        else:
            mol = Chem.MolFromSmiles(entrada)
        if mol is None:
            return jsonify({"error":"no se pudo parsear la molécula"}), 400

        # Canoniza SMILES
        smiles_canon = Chem.MolToSmiles(mol, canonical=True)

        # Opcional: SDF 3D
        mol3d = Chem.AddHs(Chem.Mol(mol))
        AllChem.EmbedMolecule(mol3d, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol3d)
        sdf = Chem.MolToMolBlock(mol3d)

        return jsonify({"smiles_canonico": smiles_canon, "sdf": sdf})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

# -------- NUEVOS ENDPOINTS QUE NECESITA EL CONSTRUCTOR 2D --------
@app.route("/smiles_to_graph", methods=["POST"])
def smiles_to_graph():
    """
    Input:  {"smiles":"c1ccccc1"}
    Output: {"atoms":[{"id":0,"el":"C","x":..,"y":..},...],
             "bonds":[{"a1":0,"a2":1,"order":1},...]}
    Notas:
    - Se generan coordenadas 2D con rdDepictor.
    - Se kekuliza para que enlaces aromáticos salgan como alternados 1/2 (útil para docencia).
    """
    data = request.get_json(silent=True) or {}
    smiles = (data.get("smiles") or "").strip()
    if not smiles:
        return jsonify({"error":"SMILES requerido"}), 400
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error":"SMILES inválido"}), 400

        # 2D coords + kekulización para evitar tipo 'AROMATIC'
        rdDepictor.Compute2DCoords(mol)
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            # si no se puede kekulizar, seguimos con lo que haya
            pass

        # coords 2D
        conf = mol.GetConformer()
        atoms = []
        for a in mol.GetAtoms():
            idx = a.GetIdx()
            pos = conf.GetAtomPosition(idx)
            atoms.append({
                "id": idx, "el": a.GetSymbol(),
                "x": float(pos.x), "y": float(pos.y)
            })

        # enlaces
        order_map = {
            Chem.rdchem.BondType.SINGLE: 1,
            Chem.rdchem.BondType.DOUBLE: 2,
            Chem.rdchem.BondType.TRIPLE: 3
        }
        bonds = []
        for b in mol.GetBonds():
            bt = b.GetBondType()
            order = order_map.get(bt, 1)  # si queda aromático u otro, lo tratamos como simple
            bonds.append({
                "a1": b.GetBeginAtomIdx(),
                "a2": b.GetEndAtomIdx(),
                "order": order
            })

        return jsonify({"atoms": atoms, "bonds": bonds})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/build_smiles", methods=["POST"])
def build_smiles():
    """
    Input:  {"atoms":[{"id":1,"el":"C",...},...],
             "bonds":[{"a1":1,"a2":2,"order":1},...]}
    Output: {"smiles":"...", "smiles_canonico":"..."}
    """
    data = request.get_json(silent=True) or {}
    atoms = data.get("atoms") or []
    bonds = data.get("bonds") or []
    if not atoms:
        return jsonify({"error":"atoms vacío"}), 400
    try:
        rw = Chem.RWMol()
        idmap = {}  # id externo -> idx RDKit
        for a in atoms:
            sym = a.get("el") or "C"
            atom = Chem.Atom(sym)
            idx = rw.AddAtom(atom)
            idmap[a.get("id")] = idx

        bond_map = {1: Chem.rdchem.BondType.SINGLE,
                    2: Chem.rdchem.BondType.DOUBLE,
                    3: Chem.rdchem.BondType.TRIPLE}
        for b in bonds:
            a1 = idmap.get(b.get("a1"))
            a2 = idmap.get(b.get("a2"))
            if a1 is None or a2 is None: 
                continue
            order = bond_map.get(int(b.get("order", 1)), Chem.rdchem.BondType.SINGLE)
            # evita duplicados (RDKit lanza si ya existe)
            if rw.GetBondBetweenAtoms(a1, a2) is None:
                rw.AddBond(a1, a2, order)

        mol = rw.GetMol()
        Chem.SanitizeMol(mol)

        smiles = Chem.MolToSmiles(mol, canonical=False)
        smiles_canon = Chem.MolToSmiles(mol, canonical=True)
        return jsonify({"smiles": smiles, "smiles_canonico": smiles_canon})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
def _http_json(url: str, timeout: int = 15):
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            if resp.status != 200:
                return None
            data = resp.read()
            return json.loads(data.decode('utf-8'))
    except Exception:
        return None

@app.route("/name_to_smiles", methods=["POST"])
def name_to_smiles():
    """
    Input:  {"query": "benzene"}  |  {"query":"cid 5793"}  |  {"query":"c1ccccc1"}
    Output: {"smiles": "<canonical smiles o el mismo texto si ya lo es>"}
    """
    data = request.get_json(silent=True) or {}
    raw = (data.get("query") or "").strip()
    if not raw:
        return jsonify({"smiles": None})

    # Si ya parece SMILES, devolver tal cual
    if re.search(r"[0-9=\[\]#()@+\-\\/]", raw):
        return jsonify({"smiles": raw})

    # ¿CID?
    m = re.match(r"^\s*(?:cid\s*[:#]?\s*)?(\d+)\s*$", raw, re.I)
    if m:
        j = _http_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{m.group(1)}/property/CanonicalSMILES/JSON")
        s = (j or {}).get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
        return jsonify({"smiles": s or raw})

    # Nombre → property
    q = urllib.parse.quote(raw)
    jA = _http_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/CanonicalSMILES/JSON")
    sA = (jA or {}).get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
    if sA:
        return jsonify({"smiles": sA})

    # Nombre → cids → property
    jC = _http_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/cids/JSON")
    cid_list = (jC or {}).get("IdentifierList", {}).get("CID", [])
    if cid_list:
        cid = cid_list[0]
        jB = _http_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON")
        sB = (jB or {}).get("PropertyTable", {}).get("Properties", [{}])[0].get("CanonicalSMILES")
        if sB:
            return jsonify({"smiles": sB})

    # Fallback
    return jsonify({"smiles": raw})

if __name__ == "__main__":
    app.run(debug=True)