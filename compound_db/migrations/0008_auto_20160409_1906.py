# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-04-09 19:06
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compound_db', '0007_auto_20160409_1902'),
    ]

    operations = [
        migrations.AlterField(
            model_name='chemblbioassaydata',
            name='operator',
            field=models.CharField(blank=True, max_length=32),
        ),
    ]
