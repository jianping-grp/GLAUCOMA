# -*- coding: utf-8 -*-
# Generated by Django 1.10.4 on 2017-04-26 09:10
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('related_info', '0006_auto_20170426_0833'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compound',
            name='drug_groups',
        ),
    ]