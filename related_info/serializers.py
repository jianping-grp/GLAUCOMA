from . import models
from dynamic_rest import serializers


class CompoundSerializer(serializers.DynamicModelSerializer):
    uniprotinfo_set = serializers.DynamicRelationField('UniprotInfoSerializer', many=True, deferred=True, embed=True)
    class Meta:
        model = models.Compound
        exclude = ['bfp', 'mol_block', 'mol']


class SynonymsSerializer(serializers.DynamicModelSerializer):
    class Meta:
        model = models.Synonyms
        exclude = []


class ProductsSerializer(serializers.DynamicModelSerializer):
    compound = serializers.DynamicRelationField('CompoundSerializer', deferred=True, embed=True)

    class Meta:
        model = models.Products
        exclude = []


class UniprotInfoSerializer(serializers.DynamicModelSerializer):
    compounds = serializers.DynamicRelationField('CompoundSerializer', many=True, deferred=True, embed=True)

    class Meta:
        model = models.UniprotInfo
        exclude = []

class KeggProteinSerializer(serializers.DynamicModelSerializer):
    keggname=serializers.DynamicRelationField('UniprotInfoSerializer')

    class Meta:
        model = models.KeggProtein
        exclude = []

class UniprotDBCompoundSerializer(serializers.DynamicModelSerializer):
    uniprot_name = serializers.DynamicRelationField('UniprotInfoSerializer')

    class Meta:
        model = models.UniprotDBCompound
        exclude = []

class UniprotAllPathway(serializers.DynamicModelSerializer):
    uniprot_name = serializers.DynamicRelationField('UniprotInfoSerializer')

    class Meta:
        model = models.UniprotAllPathway
        exclude = []
        
