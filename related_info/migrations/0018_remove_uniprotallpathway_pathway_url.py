# -*- coding: utf-8 -*-
# Generated by Django 1.11.4 on 2017-10-16 07:04
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0017_uniprotallpathway_pathway_url'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='uniprotallpathway',
            name='pathway_url',
        ),
    ]
