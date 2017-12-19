from __future__ import unicode_literals
from django.db import models
from django.utils.encoding import python_2_unicode_compatible
from django_rdkit.models import *
from django_rdkit.models.fields import *
from rdkit import Chem
from django_rdkit.config import config
from rdkit.Chem import Descriptors
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors, NumRotatableBonds

class CompoundManager(models.Manager):

    def structure_search(self, smiles, similarity):
        search_mfp2 = MORGANBV_FP(Value(smiles))
        config.tanimoto_threshold = similarity
        queryset = super(CompoundManager, self).get_queryset().filter(bfp__tanimoto=search_mfp2)
        queryset = queryset.annotate(similarity=TANIMOTO_SML('bfp', search_mfp2))
        queryset = queryset.order_by('-similarity')
        return queryset

@python_2_unicode_compatible
class Compound(models.Model):

    objects = CompoundManager()

    generic_name = models.CharField(max_length=200, blank=True)
    IUPAC_name = models.CharField(max_length=2048, blank=True)
    smiles = models.CharField(max_length=2048, blank=True)
    image = models.ImageField(upload_to='mol_images', blank=True, null=True)
    mol_weight = models.FloatField(null=True, blank=True)
    cas = models.CharField(max_length=200, blank=True, null=True)
    cid = models.CharField(max_length=200, blank=True, null=True)
    drugbank_id = models.CharField(max_length=200, blank=True, null=True)
    #drug_groups = models.CharField(max_length=200, blank=True, null=True)
    bfp = BfpField(blank=True, null=True)
    first_created = models.DateTimeField(auto_now_add=True)
    last_created = models.DateTimeField(auto_now_add=True)
    related_compounds = models.ManyToManyField('self', blank=True)
    drug_status = models.CharField(max_length=1024, blank=True, null=True)
    bindingdb_link = models.URLField(blank=True)
    chembl_link = models.URLField(blank=True)
    pathway_link = models.URLField(blank=True)

    mol = MolField(null=True, blank=True)
    mol_block = models.TextField(blank=True)
    formula = models.CharField(max_length=1024, blank=True)
    alogp = models.DecimalField(
        max_digits=9,
        decimal_places=2,
        verbose_name=('Calculated AlogP'),
        blank=True,
        null=True
    )
    hba = models.SmallIntegerField(
        verbose_name=('Number hydrogen bond acceptors'),
        blank=True,
        null=True
    )
    hbd = models.SmallIntegerField(
        verbose_name=('Number hydrogen bond donors'),
        blank=True,
        null=True
    )
    psa = models.DecimalField(
        max_digits=9,
        decimal_places=2,
        verbose_name=('Polar surface area'),
        blank=True,
        null=True
    )
    rtb = models.SmallIntegerField(
        verbose_name=('Number rotabable bonds'),
        blank=True,
        null=True
    )

    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None, *args, **kwargs):
        smiles = self.smiles
        if smiles:
            try:
                self.mol = Chem.MolFromSmiles(smiles)
                self.mol_block = Chem.MolToMolBlock(self.mol)
                self.mol_weight = Descriptors.ExactMolWt(self.mol)
                self.alogp = MolLogP(self.mol)
                self.hba = NumHAcceptors(self.mol)
                self.hbd = NumHDonors(self.mol)
                self.psa = Chem.MolSurf.TPSA(self.mol)
                self.rtb = NumRotatableBonds(self.mol)
                super(Compound, self).save(*args, **kwargs)
                self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
                self.bfp = MORGANBV_FP(Value(smiles))
            except (ValueError, TypeError):
                print "Error when storing mol object"
                pass
        super(Compound, self).save(*args, **kwargs)

    def __str__(self):
        return '{}'.format(self.generic_name or self.smiles or self.pk)

@python_2_unicode_compatible
class Synonyms(models.Model):
    name = models.CharField(max_length=1024, blank=True, null=True)
    compound = models.ForeignKey('Compound', null=True, blank=True)

    def __str__(self):
        return self.name

@python_2_unicode_compatible
class Products(models.Model):
    name = models.CharField(max_length=1024, blank=True, null=True)
    compound = models.ForeignKey('Compound', null=True, blank=True)

    def __str__(self):
        return self.name

# @python_2_unicode_compatible
# class UniprotEntry(models.Model):
#     name = models.CharField(max_length=200, blank=True)
#     compounds = models.ManyToManyField('Compound', blank=True)
#
#     def __str__(self):
#         return self.name

@python_2_unicode_compatible
class UniprotInfo(models.Model):
    entry = models.CharField(max_length=200, blank=True, null=True)
    entryname = models.CharField(max_length=200, blank=True, null=True)
    compounds = models.ManyToManyField('Compound', blank=True)
    uniprot_type = models.CharField(max_length=1024, blank=True, null=True)
    uniprot_descriptor = models.CharField(max_length=2048, blank=True, null=True)
    kegg_name = models.CharField(max_length=20, blank=True, null=True)
    kegg_url = models.URLField(max_length=1024, blank=True, null=True)
    uniprot_chembl_id = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return self.entry

    def save(self):
        self.kegg_url = 'http://www.genome.jp/dbget-bin/www_bget?{}'.format(self.kegg_name)
        super(UniprotInfo, self).save()



@python_2_unicode_compatible
class KeggProtein(models.Model):
    pdbid = models.CharField(max_length=50, blank=True, null=True)
    keggname = models.ForeignKey(UniprotInfo, blank=True, null=True)

    def __str__(self):
        return self.pdbid

@python_2_unicode_compatible
class UniprotDBCompound(models.Model):
    compound_name = models.CharField(max_length=80, blank=True, null=True)
    uniprot_name = models.ForeignKey(UniprotInfo, blank=True, null=True)
    url = models.URLField(max_length=1024, blank=True)

    def __str__(self):
        return self.compound_name

    def save(self):
        self.url = 'https://www.drugbank.ca/drugs/{}'.format(self.compound_name)
        super(UniprotDBCompound, self).save()

@python_2_unicode_compatible
class UniprotAllPathway(models.Model):
    pathway = models.CharField(max_length=1024, blank=True, null=True)
    uniprot_name = models.ForeignKey(UniprotInfo, blank=True, null=True)
    pathway_name = models.CharField(max_length=100, blank=True, null=True)
    pathway_type = models.CharField(max_length=1024, blank=True, null=True)
    pathway_url = models.URLField(max_length=1024, blank=True)

    def __str__(self):
        return self.pathway_name

    def save(self, force_insert=False, force_update=False, using=None, update_field=None):
        self.pathway_url = 'http://www.genome.jp/dbget-bin/show_pathway?{}'.format(self.pathway_name)
        super(UniprotAllPathway, self).save()


class Feedback(models.Model):
    username = models.CharField(max_length=256)
    email = models.EmailField()
    phone = models.CharField(blank=True, null=True, max_length=64)
    message = models.TextField()
    ip = models.GenericIPAddressField(blank=True, null=True)
    create_at = models.DateTimeField(auto_now=True)



