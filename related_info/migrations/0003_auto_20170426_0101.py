# -*- coding: utf-8 -*-
# Generated by Django 1.10.4 on 2017-04-26 01:01
from __future__ import unicode_literals

from django.db import migrations, models
import django_rdkit.models.fields


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0002_auto_20170425_1247'),
    ]

    operations = [
        migrations.AddField(
            model_name='compound',
            name='alogp',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=9, null=True, verbose_name='Calculated AlogP'),
        ),
        migrations.AddField(
            model_name='compound',
            name='formula',
            field=models.CharField(blank=True, max_length=1024),
        ),
        migrations.AddField(
            model_name='compound',
            name='hba',
            field=models.SmallIntegerField(blank=True, null=True, verbose_name='Number hydrogen bond acceptors'),
        ),
        migrations.AddField(
            model_name='compound',
            name='hbd',
            field=models.SmallIntegerField(blank=True, null=True, verbose_name='Number hydrogen bond donors'),
        ),
        migrations.AddField(
            model_name='compound',
            name='mol',
            field=django_rdkit.models.fields.MolField(null=True),
        ),
        migrations.AddField(
            model_name='compound',
            name='mol_block',
            field=models.TextField(blank=True),
        ),
        migrations.AddField(
            model_name='compound',
            name='psa',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=9, null=True, verbose_name='Polar surface area'),
        ),
        migrations.AddField(
            model_name='compound',
            name='rtb',
            field=models.SmallIntegerField(blank=True, null=True, verbose_name='Number rotabable bonds'),
        ),
    ]
