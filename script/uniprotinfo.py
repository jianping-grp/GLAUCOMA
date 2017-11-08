import os
import xlrd

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'GLAUCOMA.settings')
import django
django.setup()
from related_info.models import *

def upload_uniprot_info(row_number):
    print row_number
    row = table.row_values(row_number)
    entryname = row[0].strip()
    entry = row[1].strip()
    uniprot_type = row[3].strip()
    uniprot_descriptor = row[4].strip()
    kegg_name = row[5].strip()
    uniprot_chembl_id = row[7].strip()
    kegg_target = row[13].split('\n') if row[13] else []
    uniprot_compound = row[11].split('\n') if row[11] else []
    uniprot_all_pathway = row[12].split('\n') if row[12] else []


    uniprot, created = UniprotInfo.objects.get_or_create(
        entry=entry,
        entryname=entryname,
        kegg_name=kegg_name,
        uniprot_chembl_id=uniprot_chembl_id,
        uniprot_descriptor=uniprot_descriptor,
        uniprot_type=uniprot_type
    )

    for kegg_tar in kegg_target:
        keggprotein, created = KeggProtein.objects.get_or_create(pdbid=kegg_tar)
        keggprotein.keggname = uniprot
        keggprotein.save()

    for uniprotdbcom in uniprot_compound:
        uniprotcom, created = UniprotDBCompound.objects.get_or_create(compound_name=uniprotdbcom)
        uniprotcom.uniprot_name=uniprot
        uniprotcom.save()

    for uniprotpathway in uniprot_all_pathway:
        unipathway, created = UniprotAllPathway.objects.get_or_create(pathway=uniprotpathway)
        unipathway.uniprot_name = uniprot
        unipathway.save()

    # uniprot.uniprot_type = uniprot_type
    # uniprot.uniprot_descriptor = uniprot_descriptor
    # uniprot.kegg_name = kegg_name
    # uniprot.uniprot_chembl_id = uniprot_chembl_id
    # uniprot.save()



if __name__ == '__main__':
    uniprot_file = '/home/jianping/data/targets_compounds_20161111/table/new_weiyu.xlsx'
    table = xlrd.open_workbook(uniprot_file).sheet_by_index(0)
    nrows = table.nrows
    map(upload_uniprot_info, range(1, nrows))
    print 'Done!!!'
