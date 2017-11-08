from django.shortcuts import render
from django.shortcuts import HttpResponse
from collections import defaultdict
from django.core.paginator import PageNotAnInteger, EmptyPage, Paginator
from django.http import StreamingHttpResponse
from django.core import serializers
from .models import *
from django.views.generic import TemplateView, DeleteView, View
from django_rdkit.config import config
from django_rdkit.models import *
# Create your views here.

class SearchView(TemplateView):
    template_name = 'search.html'

class StructrureSearchView(View):

    def post(self, request):
        if request.is_ajax():
            data = request.POST
            is_sub = data.get("is_sub")
            tanimoto = float(data.get("tanimoto", 0.8))
            smiles = data.get("smiles")
            if is_sub and smiles:
                compounds = self._substructure_search(smiles)
                response = render(
                    request,
                    template_name='compound_result.html',
                    context={"compounds": compounds}
                )
                return StreamingHttpResponse(response.content)

            elif not is_sub and smiles:
                compounds = self._similarity_search(smiles, tanimoto=tanimoto)
                response = render(
                    request,
                    template_name='result.html',
                    context={'compounds': compounds}
                )
                return StreamingHttpResponse(response.content)


    @staticmethod
    def _similarity_search(smiles, tanimoto=0.8):
        config.tanimoto_threshold = tanimoto
        value = MORGANBV_FP(Value(smiles))
        compound_list = Compound.objects.filter(
            bfp__tanimoto=value
        )
        return compound_list


    @staticmethod
    def _substructure_search(smiles):
        compound_list = Compound.objects.filter(mol__hassubstruct=QMOL(Value(smiles)))
        return compound_list

class IdentifySearchView(View):

    def post(self, request):
        if request.is_ajax():
            data = request.POST
            type = data.get('type').strip().lower()
            query = data.get('query').strip().lower()
            context = defaultdict()
            if type == 'compound':
                compounds = self._compound_search(query)
                context.setdefault('compounds', compounds)
            elif type == 'cid':
                compounds = self._cid_search(query)
                context.setdefault('compounds', compounds)
            elif type == 'cas':
                compounds = self._cas_search(query)
                context.setdefault('compounds', compounds)
            elif type == 'formula':
                compounds = self._formula_search(query)
                context.setdefault('compounds', compounds)
            elif type == 'target':
                targets = self._target_search(query)
                context.setdefault('targets', targets)
            elif type == 'drugbank':
                compounds = self._drugbankid_search(query)
                context.setdefault('compounds', compounds)
            elif type == 'all':
                compounds = self._compound_search(query)
                targets = self._target_search(query)
                context.setdefault('compounds', compounds)
                context.setdefault('targets', targets)
            return StreamingHttpResponse(render(request, template_name='result/result.html', context=context))


    @staticmethod
    def _compound_search(query):
        compound_list = Compound.objects.filter(
            Q(generic_name__icontains=query) |
            Q(IUPAC_name__iexact=query) |
            Q(formula__iexact=query)
        )
        return compound_list

    @staticmethod
    def _target_search(query):
        target_list = UniprotInfo.object.filter(
            Q(entry__iexact=query) |
            Q(entryname__iexact=query)
        )
        return target_list

    @staticmethod
    def _formula_search(query):
        compound_list = Compound.objects.filter(
            formula__iexact=query
        )
        return compound_list

    @staticmethod
    def _cid_search(query):
        try:
            query = int(query)
            compound_list = Compound.objects.filter(
                cid__cid=query
            ).distinct()
        except ValueError:
            compound_list = Compound.objects.none()
        return compound_list

    @staticmethod
    def _cas_search(query):
        try:
            query = int(query)
            compound_list = Compound.objects.filter(
                cas__cas=query
            ).distinct()
        except ValueError:
            compound_list = Compound.objects.none()
        return compound_list

    @staticmethod
    def _drugbankid_search(query):
        try:
            query = int(query)
            compound_list = Compound.objects.filter(
                drugbank_id__drugbank_id=query
            ).distinct()
        except ValueError:
            compound_list = Compound.objects.none()
        return compound_list



class CompoundDetailView(DeleteView):
    model = Compound
    template_name = 'compound_detail.html'

class CompoundRelatedCompondsListView(TemplateView):
    template_name = 'compound_list.html'

    def get_context_data(self, **kwargs):
        context = super(CompoundRelatedCompondsListView, self).get_context_data()
        pk = int(kwargs['pk'])
        compound = Compound.objects.get(pk=pk)
        related_compounds = compound.related_compounds.all()
        context['compounds'] = related_compounds
        return context

class CompoundRelatedTargetsListView(TemplateView):
    template_name = 'target_list.html'

    def get_context_data(self, **kwargs):
        context = super(CompoundRelatedTargetsListView, self).get_context_data()
        pk = int(kwargs['pk'])
        compound = Compound.objects.get(pk=pk)
        related_targets = compound.uniprotinfo_set.all()
        context['targets'] = related_targets
        return context

class TargetDetailView(DeleteView):
    template_name = 'target_detail.html'
    model = UniprotInfo

class TargetRelatedCompoundsListView(TemplateView):
    template_name = 'compound_list.html'

    def get_context_data(self, **kwargs):
        context = super(TargetRelatedCompoundsListView, self).get_context_data()
        pk = kwargs['pk']
        target = UniprotInfo.objects.get(pk=pk)
        related_compounds = UniprotInfo.compounds.all()
        context['compounds'] = related_compounds
        return context

