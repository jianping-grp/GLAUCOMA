from django.conf.urls import url
from rest_framework import routers
from . import rest_views


rest_routers = routers.DefaultRouter()
rest_routers.register('compounds', rest_views.CompoundViewSet)
rest_routers.register('synonyms', rest_views.SynonymsViewSet)
rest_routers.register('products', rest_views.ProductsViewSet)
rest_routers.register('uniprot-info', rest_views.UniprotInfoViewSet)
rest_routers.register('kegg-protein', rest_views.KeggProteinViewSet)
rest_routers.register('uniprot-compound', rest_views.UniprotDBCompoundViewSet)
rest_routers.register('uniprot-pathway', rest_views.UniprotAllPathway)



''' from django.conf.urls import url, include
from .views import *

compounds_urls = [
    url(r'^(?P<pk>\d+)/detail/$', CompoundDetailView.as_view(), name='compound_detail'),
    url(r'^(?P<pk>\d+)/related-compounds$', CompoundRelatedCompondsListView.as_view(), name='compound_related_compounds'),
    url(r'^(?P<pk>\d+)/related-targets$', CompoundRelatedTargetsListView.as_view(), name='compund_related_targets'),
]

targets_urls = [
    url(r'^(?P<pk>\d+)/detail/$', TargetDetailView.as_view(), name='target_detail'),
    url(r'^(?P<pk>\d+)/related-compunds$', TargetRelatedCompoundsListView, name='target_related_compounds'),
]

search_urls = [
    url(r'^$', SearchView.as_view(), name='search'),
    url(r'^structure/$', StructrureSearchView.as_view(), name='structure_search'),
    url(r'^identity/$', IdentifySearchView.as_view(), name='identity_search'),
]

urlpatterns = [
    url(r'^search/', include(search_urls)),
    url(r'compounds/', include(compounds_urls)),
    url(r'targets/', include(targets_urls)),
] '''
