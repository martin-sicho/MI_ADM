from django.conf.urls import url

from compound_db.views import home, add_mol

app_name = 'compound_db'
urlpatterns = [
    url(r'^$', home, name='home')
    , url(r'^add-mol$', add_mol, name='add_mol')
]

