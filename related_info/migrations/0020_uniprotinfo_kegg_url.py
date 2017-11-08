# -*- coding: utf-8 -*-
# Generated by Django 1.11.4 on 2017-10-16 07:10
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0019_uniprotallpathway_pathway_url'),
    ]

    operations = [
        migrations.AddField(
            model_name='uniprotinfo',
            name='kegg_url',
            field=models.URLField(blank=True, max_length=1024, null=True),
        ),
    ]
