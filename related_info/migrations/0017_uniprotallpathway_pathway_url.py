# -*- coding: utf-8 -*-
# Generated by Django 1.11.4 on 2017-10-16 06:51
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0016_auto_20171016_0629'),
    ]

    operations = [
        migrations.AddField(
            model_name='uniprotallpathway',
            name='pathway_url',
            field=models.URLField(blank=True, max_length=1024),
        ),
    ]
