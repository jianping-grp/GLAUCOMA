from rest_framework import permissions, generics
from rest_framework.decorators import api_view, permission_classes, detail_route, list_route
from rest_framework import mixins
from . import models, serializers
from dynamic_rest import viewsets
from rest_framework.response import Response
from django_rdkit.models import *
from rest_framework.permissions import AllowAny
import sea

TARGET_LIST = ['CHEMBL2034',
               'CHEMBL2035',
               'CHEMBL211',
               'CHEMBL4296',
               'CHEMBL275',
               'CHEMBL1833',
               'CHEMBL1836',
               'CHEMBL3729',
               'CHEMBL3535',
               'CHEMBL3594',
               'CHEMBL3746',
               'CHEMBL3119',
               'CHEMBL261',
               'CHEMBL3912',
               'CHEMBL4619',
               'CHEMBL1881',
               'CHEMBL4884',
               'CHEMBL252',
               'CHEMBL251',
               'CHEMBL256',
               'CHEMBL254',
               'CHEMBL1951',
               'CHEMBL4235',
               'CHEMBL3510',
               'CHEMBL4425',
               'CHEMBL245',
               'CHEMBL246',
               'CHEMBL241',
               'CHEMBL1940',
               'CHEMBL1942',
               'CHEMBL3977',
               'CHEMBL230',
               'CHEMBL2056',
               'CHEMBL3710',
               'CHEMBL1867',
               'CHEMBL3242',
               'CHEMBL3969',
               'CHEMBL4789',
               'CHEMBL232',
               'CHEMBL2973',
               'CHEMBL3805',
               'CHEMBL2609',
               'CHEMBL4409',
               'CHEMBL4408',
               'CHEMBL3012',
               'CHEMBL205',
               'CHEMBL229',
               'CHEMBL222',
               'CHEMBL223',
               'CHEMBL220',
               'CHEMBL221',
               'CHEMBL226',
               'CHEMBL224',
               'CHEMBL225',
               'CHEMBL2072',
               'CHEMBL1783',
               'CHEMBL4716',
               'CHEMBL1785',
               'CHEMBL2652',
               'CHEMBL291',
               'CHEMBL290',
               'CHEMBL1916',
               'CHEMBL3025',
               'CHEMBL1914',
               'CHEMBL3878',
               'CHEMBL2885',
               'CHEMBL217',
               'CHEMBL216',
               'CHEMBL214',
               'CHEMBL213',
               'CHEMBL4640',
               'CHEMBL5932',
               'CHEMBL210',
               'CHEMBL288',
               'CHEMBL4187',
               'CHEMBL286',
               'CHEMBL3231',
               'CHEMBL2326',
               'CHEMBL1900',
               'CHEMBL1821',
               'CHEMBL5163',
               'CHEMBL1827',
               'CHEMBL5451',
               'CHEMBL3421',
               'CHEMBL1987',
               'CHEMBL4267',
               'CHEMBL1980',
               'CHEMBL2717']


class CompoundViewSet(viewsets.DynamicModelViewSet):
    queryset = models.Compound.objects.all()
    serializer_class = serializers.CompoundSerializer

    # @detail_route(methods=['POST', 'GET'])
    @list_route(methods=['POST', 'GET'], permission_classes=[permissions.AllowAny])
    def search(self, request):
        smiles = str(request.data['smiles'])
        similarity = float(request.data['similarity'])
        substructure_search = int(request.data['substructure_search'])
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


@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def target_pred(request):
    smiles = str(request.data['smiles'])
    # print smiles
    # smiles = 'CC1CCCN(C1C)C(=O)c2csc(Nc3ccc(C)cc3)n2'

    # target_list = [x['uniprot_chembl_id']
    #                for x in
    #                models.UniprotInfo.objects.filter(uniprot_chembl_id__startswith='CHEMBL').values('uniprot_chembl_id')
    #                ]
    pred_data = sea.pred2(smiles, TARGET_LIST)
    print pred_data
    return Response(pred_data)


class FeedbackCreateView(generics.CreateAPIView):
    queryset = models.Feedback.objects.all()
    serializer_class = serializers.FeedbackSerializer
    permission_classes = (AllowAny, )