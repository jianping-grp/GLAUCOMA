from rest_framework import permissions
from rest_framework.decorators import api_view, permission_classes, detail_route, list_route
from rest_framework.viewsets import ModelViewSet
from rest_framework import mixins
from . import models, serializers
from dynamic_rest import viewsets
from rest_framework.response import Response
from django_rdkit.models import *


class CompoundViewSet(viewsets.DynamicModelViewSet):
    queryset = models.Compound.objects.all()
    serializer_class = serializers.CompoundSerializer

    #@detail_route(methods=['POST', 'GET'])
    @list_route(methods=['POST', 'GET'], permission_classes=[permissions.AllowAny])
    def search(self, request):
        smiles = str(request.POST['smiles'])
        similarity = float(request.POST['similarity'])
        substructure_search = int(request.POST['substructure_search'])
        # perform substructure
        print smiles, similarity, substructure_search
        result = {}
        if substructure_search == 1:
            result = models.Compound.objects.filter(mol__hassubstruct=QMOL(Value(smiles))).all()
        # structure search
        else:
            try:
                result = models.Compound.objects.structure_search(smiles, similarity)
            except:
                print 'structure search error'
        if result:
            page = self.paginate_queryset(result)
            if page is not None:
                serializer = self.get_serializer(page, many=True)
                return self.get_paginated_response(serializer.data)
            serializer = self.get_serializer(result, many=True)
            return Response(serializer.data)

        return Response(result)


class SynonymsViewSet(viewsets.DynamicModelViewSet):
    queryset = models.Synonyms.objects.all()
    serializer_class = serializers.SynonymsSerializer


class ProductsViewSet(viewsets.DynamicModelViewSet):
    queryset = models.Products.objects.all()
    serializer_class = serializers.ProductsSerializer


class UniprotInfoViewSet(viewsets.DynamicModelViewSet):
    queryset = models.UniprotInfo.objects.all()
    serializer_class = serializers.UniprotInfoSerializer

class KeggProteinViewSet(viewsets.DynamicModelViewSet):
    queryset = models.KeggProtein.objects.all()
    serializer_class = serializers.KeggProteinSerializer

class UniprotDBCompoundViewSet(viewsets.DynamicModelViewSet):
    queryset = models.UniprotDBCompound.objects.all()
    serializer_class = serializers.UniprotDBCompoundSerializer

class UniprotAllPathway(viewsets.DynamicModelViewSet):
    queryset = models.UniprotAllPathway.objects.all()
    serializer_class = serializers.UniprotAllPathway


# @api_view(['POST'])
# @permission_classes([permissions.AllowAny])
# def structure_search(request):
#     smiles = request.POST['smiles']
#     similarity = request.POST['similarity']
#     substructure_search = int(request.POST['substructure_search'])
#     # perform substructure
#     print request.POST
#     result_json = {}
#     if substructure_search == 1:
#         pass
#     # structure search
#     else:
#         try:
#             result = models.Compound.objects.structure_search(smiles, similarity)
#         except:
#             print 'structure search error'
#
#     return Response(result)
