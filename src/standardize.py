from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize 
from rdkit.Chem.inchi import MolToInchiKey

from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

# Onstruct these objects once
un = rdMolStandardize.Uncharger()
te = rdMolStandardize.TautomerEnumerator()

failed_smiles = []

def standardize_mols(smi,):
    "Standardize molecules."
    if not isinstance(smi, str):
        failed_smiles.append({"smiles": smi, "error_type": f"Type {type(smi)} as smiles", "error": "Smiles not (str)"})
        return (None, None)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        failed_smiles.append({"smiles": smi, "error_type": "Invalid smiles", "error": "Could not convert smiles to Mol obj."})
        return (None, None)

    # Use the preferOrganic parameter for the fragmentation step
    # otherwise, this pipeline produces [Zn+2] as one of the clean smiles
    params = rdMolStandardize.CleanupParameters()
    params.preferOrganic = True

    try:
        # Remove Hs, disconnect metals (probably wouldnt need for a kinase inhibitor anyway)
        # reionize, assign stereochemistry, normalize functional groups
        mol = rdMolStandardize.Cleanup(mol)

        # Desalt / select the parent molecule
        mol = rdMolStandardize.FragmentParent(mol, params=params)

        # Neutralize
        mol = un.uncharge(mol)

        # Canonicalize Tautomers
        mol = te.Canonicalize(mol)

        # Convert to smi, inchikey
        clean_smi = Chem.MolToSmiles(mol)
        clean_inchikey = MolToInchiKey(mol)

    except Exception as e:
        failed_smiles.append({"smiles": smi, "error_type": type(e).__name__, "error": str(e)})
        return (None, None)

    return (clean_smi, clean_inchikey)