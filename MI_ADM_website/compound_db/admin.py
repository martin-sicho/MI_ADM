from django.contrib import admin

from MI_ADM_website.compound_db.models import Compound, ChEMBLTargetData, ChEMBLBioassayData

@admin.register(Compound)
class CompundAdmin(admin.ModelAdmin):
    list_display = (
        'unique_id'
        , 'smiles'
        , 'inchi_key'
        , 'description'
    )

@admin.register(ChEMBLTargetData)
class ChEMBLTargetDataAdmin(admin.ModelAdmin):
    list_display = (
        'unique_id'
        , 'uniprot_accession'
        , 'chembl_id'
        , 'organism'
        , 'preffered_name'
        , 'description'
    )

@admin.register(ChEMBLBioassayData)
class ChEMBLBioassayDataAdmin(admin.ModelAdmin):
    list_display = (
        'my_target_data'
        , 'target_data_link'
        , 'target_name'
        , 'organism'
        , 'compound_link'
        , 'value'
        , 'units'
        , 'operator'
        , 'bioactivity_type'
        , 'activity_comment'
        , 'assay_id'
    )
    search_fields = (
        'target_data__unique_id'
        , 'target_data__preffered_name'
        , 'compound__unique_id'
        , 'assay_id'
        , 'activity_comment'
    )

    def target_name(self, obj):
        return obj.target_data.preffered_name
    target_name.short_description = 'Target Name'

    def my_target_data(self, obj):
        return 'Edit'
    my_target_data.short_description = 'Edit'

    def organism(self, obj):
        return obj.target_data.organism
    organism.short_description = 'Organism'

