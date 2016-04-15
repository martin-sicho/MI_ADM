import hashlib
from datetime import datetime

from django.core.urlresolvers import reverse
from django.utils.html import escape
from django_rdkit import models
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


class Compound(models.Model):
    unique_id = models.CharField(max_length=13, unique=True, blank=False)
    description = models.TextField(blank=True)
    smiles = models.CharField(max_length=4096, blank=False, unique=True)
    inchi = models.CharField(max_length=4096, blank=False, unique=True)
    inchi_key = models.CharField(max_length=27, blank=False, unique=True)
    mol_weight_exact = models.FloatField(blank=False)
    heavy_atoms_count = models.IntegerField(blank=False)
    ring_count = models.IntegerField(blank=False)
    mol = models.MolField()

    torsionbv = models.BfpField(null=True)
    mfp2 = models.BfpField(null=True)
    ffp2 = models.BfpField(null=True)

    class MoleculeAlreadyExists(Exception):
        pass

    def _generate_id(self):
        number = datetime.now().timestamp()
        number=int(number * 10e6) # get seconds
        random_data = number.to_bytes(8, byteorder='big')
        return 'MI-M-' + hashlib.md5(random_data).hexdigest()[:8]


    def __init__(self, *args, **kwargs):
        if len(args) > 2:
            super(Compound, self).__init__(*args, **kwargs)
            return
        mol_as_RDmol = args[0] if len(args) > 0 else None
        if not mol_as_RDmol:
            mol_as_RDmol = kwargs['mol_as_RDmol'] if 'mol_as_RDmol' in kwargs else None
        if not mol_as_RDmol:
            raise RuntimeError("No RDMol specified")
        description = args[1] if len(args) > 1 else None
        if not description:
            description = kwargs['description'] if 'description' in kwargs else ''
        new_kwargs = dict()
        new_kwargs['unique_id'] = self._generate_id()
        new_kwargs['smiles'] = Chem.MolToSmiles(mol_as_RDmol)
        new_kwargs['inchi'] = Chem.MolToInchi(mol_as_RDmol)
        new_kwargs['inchi_key'] = Chem.InchiToInchiKey(new_kwargs['inchi'])
        new_kwargs['mol_weight_exact'] = Descriptors.ExactMolWt(mol_as_RDmol)
        new_kwargs['heavy_atoms_count'] = Lipinski.HeavyAtomCount(mol_as_RDmol)
        new_kwargs['ring_count'] = Lipinski.RingCount(mol_as_RDmol)
        new_kwargs['mol'] = mol_as_RDmol

        super(Compound, self).__init__(description=description, **new_kwargs)

    def save(self, *args, **kwargs):
        if Compound.objects.filter(inchi_key=self.inchi_key).exists():
            raise Compound.MoleculeAlreadyExists("Molecule with the same InchiKey was found. Cannot save.")
        super(Compound, self).save(*args, **kwargs)

    def __str__(self):
        return self.unique_id

class CompoundTrivialNames(models.Model):
    name = models.CharField(max_length=128, blank=False, unique=True)
    compound = models.ForeignKey(Compound)

class ChEMBLTargetData(models.Model):

    unique_id = models.CharField(max_length=13, unique=True, blank=False)
    uniprot_accession = models.CharField(max_length=6, blank=False)
    chembl_id = models.CharField(max_length=16, blank=False)
    organism = models.CharField(max_length=128, blank=False)
    preffered_name = models.CharField(max_length=128, blank=False)
    description = models.TextField(blank=True)

    def _generate_id(self):
        number = datetime.now().timestamp()
        number=int(number * 10e6) # get seconds
        random_data = number.to_bytes(8, byteorder='big')
        return 'MI-T-' + hashlib.md5(random_data).hexdigest()[:8]

    def __init__(self, *args, **kwargs):
        kwargs['unique_id'] = self._generate_id()
        super(ChEMBLTargetData, self).__init__(*args, **kwargs)

    def __str__(self):
        return self.unique_id

    class Meta:
        verbose_name = 'ChEMBL Target Data'
        verbose_name_plural = verbose_name

class ChEMBLBioassayData(models.Model):

    compound = models.ForeignKey(Compound, null=False)
    target_data = models.ForeignKey(ChEMBLTargetData, null=False)
    assay_id = models.CharField(max_length=32, blank=True)
    ingredient_cmpd_id = models.CharField(max_length=32, blank=True)
    units = models.CharField(max_length=32, blank=True)
    bioactivity_type = models.CharField(max_length=32, blank=True)
    value = models.FloatField(null=True)
    operator = models.CharField(max_length=32, blank=True)
    activity_comment = models.CharField(max_length=64, blank=True)
    target_confidence = models.IntegerField(null=True)

    class Meta:
        verbose_name = 'ChEMBL Bioassay Data'
        verbose_name_plural = verbose_name

    def compound_link(self):
      return '<a href="%s">%s</a>' % (reverse("admin:compound_db_compound_change", args=(self.compound.pk,)) , escape(self.compound))
    compound_link.allow_tags = True
    compound_link.short_description = 'Compound'

    def target_data_link(self):
      return '<a href="%s">%s</a>' % (reverse("admin:compound_db_chembltargetdata_change", args=(self.target_data.pk,)) , escape(self.target_data))
    target_data_link.allow_tags = True
    target_data_link.short_description = 'Target Data'