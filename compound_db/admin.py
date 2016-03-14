from django.contrib import admin

from compound_db.models import Compound, CompoundTrivialNames

admin.site.register(Compound)
admin.site.register(CompoundTrivialNames)

