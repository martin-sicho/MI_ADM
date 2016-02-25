from django.conf.urls import url

from compound_db.views import home

urlpatterns = [
    url(r'^$', home)
]

