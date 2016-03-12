import hashlib
from datetime import datetime
from django.db import models
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

class Compound(models.Model):
    unique_id = models.CharField(max_length=13, unique=True, blank=False)
    description = models.TextField(blank=True)
    smiles = models.CharField(max_length=2048, blank=False, unique=True)
    inchi = models.CharField(max_length=2048, blank=False, unique=True)
    inchi_key = models.CharField(max_length=27, blank=False, unique=True)
    mol_weight_exact = models.FloatField(blank=False)
    heavy_atoms_count = models.IntegerField(blank=False)
    ring_count = models.IntegerField(blank=False)

    class MoleculeAlreadyExists(Exception):
        pass

    def _generate_id(self):
        number = datetime.now().timestamp()
        number=int(number * 10e6) # get seconds
        random_data = number.to_bytes(8, byteorder='big')
        return 'MI-M-' + hashlib.md5(random_data).hexdigest()[:8]


    def __init__(self, mol_as_RDmol, description):
        kwargs = dict()
        kwargs['unique_id'] = self._generate_id()
        kwargs['smiles'] = Chem.MolToSmiles(mol_as_RDmol)
        kwargs['inchi'] = Chem.MolToInchi(mol_as_RDmol)
        kwargs['inchi_key'] = Chem.InchiToInchiKey(kwargs['inchi'])
        kwargs['mol_weight_exact'] = Descriptors.ExactMolWt(mol_as_RDmol)
        kwargs['heavy_atoms_count'] = Lipinski.HeavyAtomCount(mol_as_RDmol)
        kwargs['ring_count'] = Lipinski.RingCount(mol_as_RDmol)
        super(Compound, self).__init__(description=description, **kwargs)

    def save(self, *args, **kwargs):
        if Compound.objects.filter(inchi_key=self.inchi_key).exists():
            raise Compound.MoleculeAlreadyExists("Molecule with the same InchiKey was found. Cannot save.")
        super(Compound, self).save(*args, **kwargs)