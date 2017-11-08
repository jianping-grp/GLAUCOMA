# from django.contrib.auth.models import User as UserModel
# from . import models
# from graphene import AbstractType, Node
# #from graphene_django.filter import DjangoFilterConnectionField
# from graphene_django import DjangoObjectType
# import graphene
#
# class User(DjangoObjectType):
#     class Meta:
#         model = UserModel
#
# class Query(graphene.ObjectType):
#     users = graphene.List(User)
#
#     @graphene.resolve_only_args
#     def resolve_users(self):
#         return UserModel.objects.all()

# schema = graphene.Schema(query=Query)

# from . import models
# from graphene import AbstractType, Node
# from graphene_django.filter import DjangoFilterConnectionField
# from graphene_django.types import DjangoObjectType
#
# class CompoundNode(DjangoObjectType):
#     class Meta:
#         model = models.Compound
#         interfaces = (Node, )
#         exclude_fields = [
#             'mol',
#         ]
#         filter_fields = [
#             'generic_name',
#             'cid',
#             'cas',
#             'drugbank_id',
#         ]
#
#     @classmethod
#     def get_node(cls, id, context, info):
#         try:
#             cpd = cls._meta.model.objects.get(id=id)
#         except cls._meta.model.DoesNotExist:
#             return None
#         return cpd
#
# class SynonymsNode(DjangoObjectType):
#     class Meta:
#         model = models.Synonyms
#         interfaces = (Node, )
#
# class ProductsNode(DjangoObjectType):
#     class Meta:
#         model = models.Products
#         interfaces = (Node, )
#
# class UniprotInfoNode(DjangoObjectType):
#     class Meta:
#         model = models.UniprotInfo
#         filter_fields = ['entry', 'entryname', 'compounds']
#         interfaces = (Node, )
#
#
#
# class Query(AbstractType):
#     compound = Node.Field(CompoundNode)
#     all_compound = DjangoFilterConnectionField(CompoundNode)
#
#     Synonym = Node.Field(SynonymsNode)
#     all_synonyms = DjangoFilterConnectionField(SynonymsNode)
#
#     product = Node.Field(ProductsNode)
#     all_products = DjangoFilterConnectionField(ProductsNode)
#
#     uniprotinfo = Node.Field(UniprotInfoNode)
#     all_uniprotinfo = DjangoFilterConnectionField(UniprotInfoNode)

from .models import *
from django_filters import FilterSet, OrderingFilter
from graphene import ObjectType, Schema, relay, AbstractType
from graphene_django import DjangoObjectType
from graphene_django.filter import DjangoFilterConnectionField


class CompoundNode(DjangoObjectType):
    class Meta:
        model = Compound
        interfaces = (relay.Node, )
        exclude_fields = ['mol']
        filter_fields = {
            'generic_name': ['exact', 'icontains', 'istartswith']
        }


class UniprotInfoNode(DjangoObjectType):
    compounds = DjangoFilterConnectionField(CompoundNode)
    class Meta:
        model = UniprotInfo
        interfaces = (relay.Node, )
        #filter_fields = {
        #    'entry': ['exact'],
        #    'entryname': ['exact']
        #}


class SynonymNode(DjangoObjectType):
    class Meta:
        model = Synonyms
        interfaces = (relay.Node, )
        exclude_fields = []

class ProductNode(DjangoObjectType):
    class Meta:
        model = Products
        interfaces = (relay.Node, )
        exclude_fields = []


class Query(AbstractType):
    compound = relay.Node.Field(CompoundNode)
    all_compounds = DjangoFilterConnectionField(CompoundNode)

    uniprot = relay.Node.Field(UniprotInfoNode)
    all_uniprots = DjangoFilterConnectionField(UniprotInfoNode)

    product = relay.Node.Field(ProductNode)
    all_products = DjangoFilterConnectionField(ProductNode)

    synonym = relay.Node.Field(SynonymNode)
    all_synonyms = DjangoFilterConnectionField(SynonymNode)





