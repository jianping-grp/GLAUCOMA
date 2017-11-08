import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'GLAUCOMA.settings')
import django
django.setup()

from related_info.models import *

# def upload_compound_chembl_info(index):
#     compound = Compound.objects.get(pk=index)
#     if compound.generic_name.startswith()

compound_list = list(Compound.objects.filter(generic_name__icontains='CHEMBL'))
for compound in compound_list:
    print compound
    chembl_url = 'https://www.ebi.ac.uk/chembl/compound/inspect/{}'.format(compound.generic_name)
    compound.chembl_link = chembl_url
    compound.save()

