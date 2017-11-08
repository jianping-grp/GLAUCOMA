import os
import xlrd

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'GLAUCOMA.settings')
import django
django.setup()
from related_info.models import *
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned
import logging

logging.basicConfig(
    level=logging.WARNING,
    format="%(levelname)s-[%(message)s]",
    filename="/home/jianping/django_test/logs/compound_log.txt",
    filemode='w'
)

def uniprot_extract(cell):
    herbs = cell.split('\n') if cell else []
    for herb in herbs:
        yield herb.split(':')

def upload_compound(row_number):
    row = table.row_values(row_number)
    smile = row[0].strip()
    related_products = row[1].split('\n') if row[1] else []
    synonyms = row[2].split('\n') if row[2] else []
    iupac_name = row[3].strip()
    drugbank_id = row[4].strip()
    drug_groups = row[5].strip()
    generic_name = row[6].strip()
    #uniprot_entrys = row[8].split('\n') if row[8] else []
    cid = row[9].strip()
    chemblid = row[10].strip()
    uniprot_info = uniprot_extract(row[12])

    compound, created = Compound.objects.get_or_create(
        smiles=smile,
        IUPAC_name=iupac_name,

        drugbank_id=drugbank_id,
        drug_status=drug_groups,
        generic_name=generic_name,
        cid=cid,
    )

    for related_product in related_products:
        try:
            product, created = Products.objects.get_or_create(name=related_product)
            product.compound = compound
            product.save()
        except Products.DoesNotExist:
            logging.warning("{} Does not find in Product database".format(related_product))
        except Products.MultipleObjectsReturned:
            logging.warning("{} find multiple objs in database".format(related_product))


    for synonym in synonyms:
        try:
            snnm, created = Synonyms.objects.get_or_create(name=synonym)
            snnm.compound = compound
            snnm.save()
        except Synonyms.DoesNotExist:
            logging.warning("{} Does not find in Synonyms database".format(synonym))
        except Synonyms.MultipleObjectsReturned:
            logging.warning("{} find multiple synonym objs in database".format(synonym))

    for unip in uniprot_info:
        print unip
        try:
            u, created = UniprotInfo.objects.get_or_create(entryname=unip[0], entry=unip[1])
            u.compounds.add(compound)
            u.save()
        except UniprotInfo.DoesNotExist:
            logging.info("Can not find %s in database" % unip)
        except UniprotInfo.MultipleObjectsReturned:
            logging.info("Return multiple uniprot_entryname")


if __name__ == '__main__':
    compound_file = '/home/jianping/data/targets_compounds_20161111/table/integrated_data_2.xlsx'
    table = xlrd.open_workbook(compound_file).sheet_by_index(0)
    nrows = table.nrows
    map(upload_compound, range(1, nrows))
    print 'Done!!!'

