import graphene
from graphene_django.debug import DjangoDebug
import related_info.schema

class Query(related_info.schema.Query, graphene.ObjectType):
    debug = graphene.Field(DjangoDebug, name='__debug')

schema = graphene.Schema(query=Query)
