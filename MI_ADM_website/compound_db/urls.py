from django.conf.urls import url

from compound_db.views import home, add_mol, autocomplete_api, search_api, molimages_api, add_chembl_mols

app_name = 'compound_db'
urlpatterns = [
    url(r'^$', home, name='home')
    , url(r'^add-mol$', add_mol, name='add_mol')
    , url(r'^add-chembl-mols$', add_chembl_mols, name='add_chembl_mols')

    # database API URLs
    , url(r'^api/search$', search_api, name='search_api')
    , url(r'^api/autocomplete$', autocomplete_api, name='autocomplete_api')
    , url(r'^api/molimages/(?P<unique_id>MI-M-.{8})', molimages_api, name='molimages_api')
]

