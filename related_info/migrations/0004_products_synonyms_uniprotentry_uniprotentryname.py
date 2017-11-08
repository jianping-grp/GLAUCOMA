# -*- coding: utf-8 -*-
# Generated by Django 1.10.4 on 2017-04-26 01:10
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0003_auto_20170426_0101'),
    ]

    operations = [
        migrations.CreateModel(
            name='Products',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=1024, null=True)),
                ('compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='related_info.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='Synonyms',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=1024, null=True)),
                ('compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='related_info.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='UniprotEntry',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=200)),
                ('compounds', models.ManyToManyField(blank=True, to='related_info.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='UniprotEntryName',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=200, null=True)),
                ('uniprotentry', models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='related_info.UniprotEntry')),
            ],
        ),
    ]
