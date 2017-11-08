# -*- coding: utf-8 -*-
# Generated by Django 1.10.4 on 2017-04-26 08:32
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0004_products_synonyms_uniprotentry_uniprotentryname'),
    ]

    operations = [
        migrations.CreateModel(
            name='UniprotInfo',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('entry', models.CharField(blank=True, max_length=200, null=True)),
                ('entryname', models.CharField(blank=True, max_length=200, null=True)),
                ('compounds', models.ManyToManyField(blank=True, null=True, to='related_info.Compound')),
            ],
        ),
        migrations.RemoveField(
            model_name='uniprotentry',
            name='compounds',
        ),
        migrations.RemoveField(
            model_name='uniprotentryname',
            name='uniprotentry',
        ),
        migrations.DeleteModel(
            name='UniprotEntry',
        ),
        migrations.DeleteModel(
            name='UniprotEntryName',
        ),
    ]
