from django.contrib import admin

from compound_db.models import Compound, CompoundTrivialNames, ChEMBLTargetData, ChEMBLBioassayData

admin.site.register(Compound)
admin.site.register(CompoundTrivialNames)
admin.site.register(ChEMBLTargetData)
admin.site.register(ChEMBLBioassayData)

