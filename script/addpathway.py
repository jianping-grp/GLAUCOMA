import os

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'GLAUCOMA.settings')
import django
django.setup()

from related_info.models import *

# print len(UniprotAllPathway.objects.all())

# pathways = UniprotAllPathway.objects.all()
# print pathways
# for pathway in pathways:
#     print pathway
    # temp_list = pathway.pathway.split('#')
    # pathway.pathway_name = temp_list[0]
    # pathway.pathway_type = temp_list[1]
    # pathway.save()

def add_url_info(index):
    compound = UniprotDBCompound.objects.get(pk=idx)
    # temp_list = pathway.pathway.split('#')
    # pathway.pathway_name = temp_list[0]
    # pathway.pathway_type = temp_list[1]
    compound.save()

el_num = len(UniprotDBCompound.objects.all())
for idx in range(1, el_num):
    print idx
    add_url_info(idx)


